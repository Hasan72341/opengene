from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from typing import Sequence


APP_NAME = "OpenGene Studio"
APP_VERSION = "0.1.0"
_DESKTOP_DIR = Path(__file__).resolve().parent
_QML_ENTRYPOINT = _DESKTOP_DIR / "qml" / "Main.qml"


def qml_entrypoint() -> Path:

    return _QML_ENTRYPOINT


def pyside6_is_available() -> bool:

    return importlib.util.find_spec("PySide6") is not None


def missing_runtime_message() -> str:

    return "\n".join(
        (
            "PySide6 is not installed, so the optional native desktop host cannot start.",
            "",
            "The touch-first browser workstation remains available:",
            "  npm run dev",
            "  # then open http://localhost:4174",
            "",
            "To install this optional host and its Qt runtime:",
            "  python -m pip install -e ./apps/desktop",
            "  opengene-desktop",
            "",
            "No scientific operation is performed by this fallback message.",
        )
    )


def _core_status() -> str:

    if importlib.util.find_spec("opengene_core") is None:
        return "UI preview — portable molecular core is not installed"
    return "Portable molecular core detected"


def _controller_type(QtCore: object) -> type:

    QObject = QtCore.QObject
    Property = QtCore.Property
    Signal = QtCore.Signal
    Slot = QtCore.Slot

    class AppController(QObject):

        statusMessageChanged = Signal()

        def __init__(self) -> None:
            super().__init__()
            self._status_message = "Native host ready — choose a workspace"
            self._core_status = _core_status()

        @Property(str, constant=True)
        def appName(self) -> str:  # noqa: N802 - Qt property names are QML API
            return APP_NAME

        @Property(str, constant=True)
        def appVersion(self) -> str:  # noqa: N802
            return APP_VERSION

        @Property(str, constant=True)
        def coreStatus(self) -> str:  # noqa: N802
            return self._core_status

        @Property(str, notify=statusMessageChanged)
        def statusMessage(self) -> str:  # noqa: N802
            return self._status_message

        @Slot(str)
        def announce(self, message: str) -> None:

            message = str(message).strip()
            if message and message != self._status_message:
                self._status_message = message
                self.statusMessageChanged.emit()

    return AppController


def launch(qt_arguments: Sequence[str] | None = None) -> int:

    if not pyside6_is_available():
        print(missing_runtime_message(), file=sys.stderr)
        return 2

    from PySide6.QtCore import QCoreApplication, QUrl  # type: ignore[import-not-found]
    from PySide6.QtGui import QGuiApplication  # type: ignore[import-not-found]
    from PySide6.QtQml import QQmlApplicationEngine  # type: ignore[import-not-found]
    from PySide6 import QtCore  # type: ignore[import-not-found]

    if not qml_entrypoint().is_file():
        print(
            f"Native host entry point is missing: {qml_entrypoint()}", file=sys.stderr
        )
        return 1

    QCoreApplication.setOrganizationName("OpenGene")
    QCoreApplication.setOrganizationDomain("opengene.local")
    QCoreApplication.setApplicationName(APP_NAME)
    QCoreApplication.setApplicationVersion(APP_VERSION)

    app = QGuiApplication([sys.argv[0], *(qt_arguments or ())])
    engine = QQmlApplicationEngine()
    controller = _controller_type(QtCore)()
    engine.rootContext().setContextProperty("appController", controller)
    engine.addImportPath(str(qml_entrypoint().parent))
    engine.load(QUrl.fromLocalFile(str(qml_entrypoint())))

    if not engine.rootObjects():
        print("OpenGene native host could not load its QML shell.", file=sys.stderr)
        return 1

    return app.exec()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Launch the optional PySide6 / Qt Quick OpenGene desktop host."
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="report whether the optional PySide6 runtime is available",
    )
    parser.add_argument(
        "--qml-path",
        action="store_true",
        help="print the installed QML entry point and exit",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {APP_VERSION}")
    return parser


def main(argv: Sequence[str] | None = None) -> int:

    parser = build_parser()
    parsed, qt_arguments = parser.parse_known_args(argv)

    if parsed.qml_path:
        print(qml_entrypoint())
        return 0

    if parsed.check:
        if pyside6_is_available():
            print("PySide6 runtime available; native OpenGene host can launch.")
            return 0
        print(missing_runtime_message(), file=sys.stderr)
        return 1

    return launch(qt_arguments)


if __name__ == "__main__":  # pragma: no cover - direct local invocation
    raise SystemExit(main())
