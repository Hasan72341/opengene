import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: tabBar

    property var colors
    property int currentIndex: 0
    signal destinationRequested(int index)

    implicitHeight: colors.compactTouchTarget + 14
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

    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 8
        anchors.rightMargin: 8
        spacing: 4

        Repeater {
            model: tabBar.destinations

            delegate: NavigationButton {
                required property var modelData
                required property int index

                Layout.fillWidth: true
                Layout.fillHeight: true
                colors: tabBar.colors
                compact: true
                label: modelData.label
                glyph: modelData.glyph
                selected: tabBar.currentIndex === index
                Accessible.name: modelData.label + " workspace"
                onClicked: tabBar.destinationRequested(index)
            }
        }
    }
}
