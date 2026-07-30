import QtQuick
import QtQuick.Controls

Button {
    id: control

    property bool selected: false
    property string glyph: ""
    property string label: ""
    property var colors
    property bool compact: false

    implicitWidth: compact ? 72 : 176
    implicitHeight: compact ? colors.compactTouchTarget : colors.touchTarget
    padding: 8
    hoverEnabled: !compact

    background: Rectangle {
        radius: 12
        color: control.down
               ? colors.accentSoft
               : (control.selected ? colors.accentSoft : "transparent")
        border.width: control.selected ? 1 : 0
        border.color: colors.accent
    }

    contentItem: Row {
        spacing: 10
        anchors.centerIn: parent

        Text {
            text: control.glyph
            font.pixelSize: 20
            color: control.selected ? colors.accent : colors.mutedInk
            horizontalAlignment: Text.AlignHCenter
            verticalAlignment: Text.AlignVCenter
        }

        Text {
            visible: !control.compact
            text: control.label
            font.family: colors.uiFont
            font.pixelSize: 14
            font.weight: control.selected ? Font.DemiBold : Font.Medium
            color: control.selected ? colors.ink : colors.mutedInk
            elide: Text.ElideRight
            verticalAlignment: Text.AlignVCenter
        }
    }
}
