from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest

from apps.cli.main import main


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def run_cli(*arguments: str) -> tuple[int, dict[str, object], str]:
    output = StringIO()
    errors = StringIO()
    with redirect_stdout(output), redirect_stderr(errors):
        code = main(list(arguments))
    return code, json.loads(output.getvalue()), errors.getvalue()


class CliWorkflowTests(unittest.TestCase):
    def test_alignment_and_pcr_are_scriptable(self) -> None:
        sequence = str(FIXTURES / "cli-demo.fasta")
        alignment_code, alignment, alignment_error = run_cli("align", sequence, sequence)
        self.assertEqual(alignment_code, 0, alignment_error)
        self.assertEqual(alignment["cigar"], "42=")

        pcr_code, pcr, pcr_error = run_cli(
            "pcr",
            sequence,
            "GAATTCATGAAACCCGGG",
            "GGATCCTTTAAACCCGGG",
        )
        self.assertEqual(pcr_code, 0, pcr_error)
        self.assertTrue(pcr["feasible"])
        self.assertEqual(len(pcr["products"]), 1)

    def test_verified_ogs_project_round_trips_through_cli(self) -> None:
        sequence = str(FIXTURES / "cli-demo.fasta")
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = str(Path(temporary_directory) / "demo.ogs")
            create_code, created, create_error = run_cli(
                "project", "create", output, "CLI test project", sequence
            )
            self.assertEqual(create_code, 0, create_error)
            self.assertTrue(Path(output).is_file())
            self.assertEqual(created["document_count"], 1)

            info_code, info, info_error = run_cli("project", "info", output)
            self.assertEqual(info_code, 0, info_error)
            self.assertEqual(info["name"], "CLI test project")
            self.assertEqual(info["documents"][0]["length"], 42)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
