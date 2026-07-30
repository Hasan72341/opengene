import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: rail

    property var colors
    property int currentIndex: 0
    signal destinationRequested(int index)

    implicitWidth: 208
    color: colors.surface
    border.width: 1
    border.color: colors.border

    readonly property var destinations: [
        { "label": "Map", "glyph": "◉" },
        { "label": "Sequence", "glyph": "≡" },
        { "label": "Features", "glyph": "⌕" },
        { "label": "Tools", "glyph": "⊞" },
        { "label": "History", "glyph": "◌" }
    ]

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 14
        spacing: 8

        Label {
            text: "WORKSPACES"
            color: rail.colors.mutedInk
            font.pixelSize: 11
            font.letterSpacing: 1.1
            Layout.leftMargin: 10
            Layout.topMargin: 6
        }

        Repeater {
            model: rail.destinations

            delegate: NavigationButton {
                required property var modelData
                required property int index

                Layout.fillWidth: true
                colors: rail.colors
                label: modelData.label
                glyph: modelData.glyph
                selected: rail.currentIndex === index
                onClicked: rail.destinationRequested(index)
            }
        }

        Item { Layout.fillHeight: true }

        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: 84
            radius: 12
            color: rail.colors.surfaceMuted

            Column {
                anchors.fill: parent
                anchors.margins: 12
                spacing: 5

                Text {
                    text: "Local-first workspace"
                    color: rail.colors.ink
                    font.pixelSize: 12
                    font.weight: Font.DemiBold
                }
                Text {
                    text: "No cloud account required"
                    color: rail.colors.mutedInk
                    font.pixelSize: 11
                    wrapMode: Text.WordWrap
                    width: parent.width
                }
            }
        }
    }
}
