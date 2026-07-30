from __future__ import annotations

import argparse
import ast
import io
import re
import sys
import tokenize
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SKIP_DIRS = {
    ".git",
    "node_modules",
    "__pycache__",
    "build",
    "dist",
    ".venv",
    "venv",
    ".pytest_cache",
    "playwright-report",
    "test-results",
    "scripts",
}

KEEP_PYTHON_COMMENT = re.compile(
    r"^#\s*(?:noqa\b|type:\s*ignore|pragma\b|pylint:|mypy:|ruff:|fmt:\s*(?:on|off)|isort:|-\*-\s*coding|coding[:=])"
)

JS_LIKE_SUFFIXES = {".js", ".mjs", ".cjs", ".qml"}
CSS_SUFFIXES = {".css"}
HTML_SUFFIXES = {".html", ".htm"}
PYTHON_SUFFIXES = {".py"}

REGEX_PRECEDING = set("(,=:[!&|?{};+-*%~^<>")
REGEX_KEYWORDS = {
    "return",
    "typeof",
    "instanceof",
    "in",
    "of",
    "new",
    "delete",
    "void",
    "do",
    "else",
    "case",
    "yield",
    "await",
    "throw",
}


COMMENT_MARK = "\x00"


@dataclass
class Result:
    path: Path
    changed: bool
    removed_lines: int
    error: str | None = None


def _finalize(text: str) -> str:
    lines: list[str] = []
    for line in text.split("\n"):
        if COMMENT_MARK in line and not line.replace(COMMENT_MARK, "").strip():
            continue
        lines.append(line.replace(COMMENT_MARK, "").rstrip())
    result = "\n".join(_collapse_blank_runs(lines))
    return result + "\n" if result else result


def _erased(comment: str) -> str:
    return COMMENT_MARK + f"\n{COMMENT_MARK}" * comment.count("\n")


def _collapse_blank_runs(lines: list[str]) -> list[str]:
    output: list[str] = []
    blanks = 0
    for line in lines:
        if line.strip():
            blanks = 0
            output.append(line)
            continue
        blanks += 1
        if blanks <= 2 and output:
            output.append(line)
    while output and not output[-1].strip():
        output.pop()
    while output and not output[0].strip():
        output.pop(0)
    return output


def _docstring_nodes(tree: ast.AST) -> list[tuple[ast.Expr, bool]]:
    found: list[tuple[ast.Expr, bool]] = []
    holders = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    for node in ast.walk(tree):
        if not isinstance(node, holders):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        first = body[0]
        if not isinstance(first, ast.Expr):
            continue
        value = first.value
        if not (isinstance(value, ast.Constant) and isinstance(value.value, str)):
            continue
        only_statement = len(body) == 1 and not isinstance(node, ast.Module)
        found.append((first, only_statement))
    return found


def _strip_docstrings_from_tree(tree: ast.AST) -> ast.AST:
    for node, _ in _docstring_nodes(tree):
        for parent in ast.walk(tree):
            body = getattr(parent, "body", None)
            if isinstance(body, list) and body and body[0] is node:
                body.pop(0)
                if not body and not isinstance(parent, ast.Module):
                    body.append(ast.Pass())
    return ast.fix_missing_locations(tree)


def strip_python(source: str) -> str:
    tree = ast.parse(source)
    lines = source.splitlines()
    drop: set[int] = set()
    replace: dict[int, str] = {}

    for node, only_statement in _docstring_nodes(tree):
        start = node.lineno - 1
        end = node.end_lineno - 1
        for index in range(start, end + 1):
            drop.add(index)
        if only_statement:
            indent = " " * (node.col_offset)
            replace[start] = f"{indent}pass"

    tokens = tokenize.generate_tokens(io.StringIO(source).readline)
    comment_positions: dict[int, list[int]] = {}
    for token in tokens:
        if token.type != tokenize.COMMENT:
            continue
        row = token.start[0] - 1
        if row in drop:
            continue
        text = token.string
        if row == 0 and text.startswith("#!"):
            continue
        if KEEP_PYTHON_COMMENT.match(text):
            continue
        comment_positions.setdefault(row, []).append(token.start[1])

    for row, columns in comment_positions.items():
        line = lines[row]
        for column in sorted(columns, reverse=True):
            line = line[:column]
        stripped = line.rstrip()
        if not stripped.strip():
            drop.add(row)
        else:
            lines[row] = stripped

    output: list[str] = []
    for index, line in enumerate(lines):
        if index in replace:
            output.append(replace[index])
            continue
        if index in drop:
            continue
        output.append(line)

    result = "\n".join(_collapse_blank_runs(output))
    if result:
        result += "\n"

    expected = ast.dump(_strip_docstrings_from_tree(ast.parse(source)))
    actual = ast.dump(_strip_docstrings_from_tree(ast.parse(result)))
    if expected != actual:
        raise ValueError("python AST changed after stripping")
    return result


def strip_js_like(source: str) -> str:
    output: list[str] = []
    index = 0
    length = len(source)
    template_stack: list[int] = []
    previous_significant = ""
    previous_word = ""

    def regex_allowed() -> bool:
        if previous_word in REGEX_KEYWORDS:
            return True
        if previous_significant == "":
            return True
        return previous_significant in REGEX_PRECEDING

    while index < length:
        char = source[index]
        pair = source[index : index + 2]

        if pair == "//":
            end = source.find("\n", index)
            output.append(COMMENT_MARK)
            index = length if end == -1 else end
            continue

        if pair == "/*":
            end = source.find("*/", index + 2)
            if end == -1:
                index = length
                continue
            output.append(_erased(source[index : end + 2]))
            index = end + 2
            continue

        if char in "'\"":
            quote = char
            output.append(char)
            index += 1
            while index < length:
                current = source[index]
                output.append(current)
                index += 1
                if current == "\\" and index < length:
                    output.append(source[index])
                    index += 1
                    continue
                if current == quote:
                    break
            previous_significant = quote
            previous_word = ""
            continue

        if char == "`":
            output.append(char)
            index += 1
            template_stack.append(0)
            while index < length and template_stack:
                current = source[index]
                if current == "\\" and index + 1 < length:
                    output.append(source[index : index + 2])
                    index += 2
                    continue
                if current == "`":
                    output.append(current)
                    index += 1
                    template_stack.pop()
                    break
                if source[index : index + 2] == "${":
                    output.append("${")
                    index += 2
                    depth = 1
                    while index < length and depth:
                        inner = source[index]
                        if inner == "{":
                            depth += 1
                        elif inner == "}":
                            depth -= 1
                            if depth == 0:
                                output.append("}")
                                index += 1
                                break
                        output.append(inner)
                        index += 1
                    continue
                output.append(current)
                index += 1
            previous_significant = "`"
            previous_word = ""
            continue

        if char == "/" and regex_allowed():
            scan = index + 1
            in_class = False
            terminated = False
            while scan < length:
                current = source[scan]
                if current == "\\":
                    scan += 2
                    continue
                if current == "\n":
                    break
                if current == "[":
                    in_class = True
                elif current == "]":
                    in_class = False
                elif current == "/" and not in_class:
                    terminated = True
                    break
                scan += 1
            if terminated:
                while scan + 1 < length and source[scan + 1].isalpha():
                    scan += 1
                output.append(source[index : scan + 1])
                index = scan + 1
                previous_significant = "/"
                previous_word = ""
                continue

        output.append(char)
        if not char.isspace():
            previous_significant = char
            if char.isalnum() or char in "_$":
                previous_word += char
            else:
                previous_word = ""
        elif char == "\n":
            previous_word = ""
        index += 1

    return _finalize("".join(output))


def strip_css(source: str) -> str:
    output: list[str] = []
    index = 0
    length = len(source)
    while index < length:
        char = source[index]
        if source[index : index + 2] == "/*":
            end = source.find("*/", index + 2)
            if end == -1:
                index = length
                continue
            output.append(_erased(source[index : end + 2]))
            index = end + 2
            continue
        if char in "'\"":
            quote = char
            output.append(char)
            index += 1
            while index < length:
                current = source[index]
                output.append(current)
                index += 1
                if current == "\\" and index < length:
                    output.append(source[index])
                    index += 1
                    continue
                if current == quote:
                    break
            continue
        output.append(char)
        index += 1

    return _finalize("".join(output))


def strip_html(source: str) -> str:
    output: list[str] = []
    index = 0
    length = len(source)
    lowered = source.lower()
    while index < length:
        if lowered.startswith("<!--", index):
            end = source.find("-->", index + 4)
            if end == -1:
                index = length
                continue
            output.append(_erased(source[index : end + 3]))
            index = end + 3
            continue
        for tag, stripper in (("script", strip_js_like), ("style", strip_css)):
            open_tag = f"<{tag}"
            if not lowered.startswith(open_tag, index):
                continue
            head_end = source.find(">", index)
            close = lowered.find(f"</{tag}", head_end + 1)
            if head_end == -1 or close == -1:
                continue
            head = source[index : head_end + 1]
            body = source[head_end + 1 : close]
            output.append(head)
            output.append(stripper(body) if body.strip() else body)
            index = close
            break
        else:
            output.append(source[index])
            index += 1
            continue

    return _finalize("".join(output))


STRIPPERS = [
    (PYTHON_SUFFIXES, strip_python),
    (JS_LIKE_SUFFIXES, strip_js_like),
    (CSS_SUFFIXES, strip_css),
    (HTML_SUFFIXES, strip_html),
]


def stripper_for(path: Path):
    for suffixes, stripper in STRIPPERS:
        if path.suffix in suffixes:
            return stripper
    return None


def iter_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if stripper_for(path):
            yield path


def process(path: Path, root: Path, write: bool) -> Result:
    stripper = stripper_for(path)
    original = path.read_text(encoding="utf-8")
    try:
        updated = stripper(original)
    except Exception as error:
        return Result(path.relative_to(root), False, 0, f"{type(error).__name__}: {error}")
    if updated == original:
        return Result(path.relative_to(root), False, 0)
    removed = len(original.splitlines()) - len(updated.splitlines())
    if write:
        path.write_text(updated, encoding="utf-8")
    return Result(path.relative_to(root), True, removed)


def main() -> int:
    parser = argparse.ArgumentParser(description="Strip comments and docstrings from the OpenGene tree")
    parser.add_argument("--root", default=str(ROOT), help="directory to process")
    parser.add_argument("--check", action="store_true", help="report without writing")
    arguments = parser.parse_args()

    root = Path(arguments.root).resolve()
    results = [process(path, root, not arguments.check) for path in iter_files(root)]
    errors = [result for result in results if result.error]
    changed = [result for result in results if result.changed]

    for result in changed:
        print(f"{'would strip' if arguments.check else 'stripped'} {result.path} (-{result.removed_lines} lines)")
    for result in errors:
        print(f"ERROR {result.path}: {result.error}", file=sys.stderr)

    total = sum(result.removed_lines for result in changed)
    print(f"\n{len(changed)}/{len(results)} files changed, {total} lines removed")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
