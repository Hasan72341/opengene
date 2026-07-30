import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: card

    property var colors
    property string title: ""
    property string detail: ""
    property string actionLabel: "Open"
    property bool actionEnabled: true
    signal actionRequested()

    implicitHeight: content.implicitHeight + 32
    radius: colors.cornerRadius
    color: colors.surface
    border.width: 1
    border.color: colors.border

    ColumnLayout {
        id: content
        anchors.fill: parent
        anchors.margins: 16
        spacing: 10

        Label {
            text: card.title
            color: card.colors.ink
            font.pixelSize: 17
            font.weight: Font.DemiBold
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
        }

        Label {
            text: card.detail
            color: card.colors.mutedInk
            font.pixelSize: 13
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
        }

        Button {
            Layout.topMargin: 2
            Layout.minimumHeight: card.colors.touchTarget
            Layout.minimumWidth: 132
            enabled: card.actionEnabled
            text: card.actionLabel
            onClicked: card.actionRequested()

            background: Rectangle {
                radius: 10
                color: parent.down
                       ? card.colors.accentSoft
                       : (parent.enabled ? card.colors.accent : card.colors.border)
            }
            contentItem: Text {
                text: parent.text
                color: parent.enabled ? "white" : card.colors.mutedInk
                horizontalAlignment: Text.AlignHCenter
                verticalAlignment: Text.AlignVCenter
                font.pixelSize: 14
                font.weight: Font.DemiBold
            }
        }
    }
}
