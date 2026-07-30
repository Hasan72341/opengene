import QtQuick
import QtQuick.Controls

Item {
    id: preview

    property var colors
    property real zoom: 1.0
    property string selectedFeature: "Circular construct"
    signal featureActivated(string label)

    implicitHeight: 440
    clip: true

    Rectangle {
        anchors.fill: parent
        radius: preview.colors.cornerRadius
        color: preview.colors.surface
        border.width: 1
        border.color: preview.colors.border
    }

    Item {
        id: plot
        width: Math.min(preview.width * 0.78, preview.height * 0.78)
        height: width
        anchors.centerIn: parent
        scale: preview.zoom

        Canvas {
            id: plasmidCanvas
            anchors.fill: parent
            antialiasing: true
            onPaint: {
                var context = getContext("2d")
                context.clearRect(0, 0, width, height)

                var centerX = width / 2
                var centerY = height / 2
                var radius = Math.min(width, height) * 0.30
                var ringWidth = Math.max(19, width * 0.06)

                context.lineWidth = ringWidth
                context.lineCap = "round"

                context.strokeStyle = preview.colors.plasmid
                context.beginPath()
                context.arc(centerX, centerY, radius, 0, Math.PI * 2)
                context.stroke()

                drawArc(context, centerX, centerY, radius, ringWidth + 7,
                        -1.45, 0.34, preview.colors.featureOne)
                drawArc(context, centerX, centerY, radius, ringWidth + 7,
                        0.55, 1.88, preview.colors.featureTwo)
                drawArc(context, centerX, centerY, radius, ringWidth + 7,
                        2.12, 3.55, preview.colors.featureThree)

                context.fillStyle = preview.colors.ink
                context.font = Math.max(13, width * 0.055) + "px Sans Serif"
                context.textAlign = "center"
                context.fillText("6,214 bp", centerX, centerY - 3)
                context.fillStyle = preview.colors.mutedInk
                context.font = Math.max(11, width * 0.04) + "px Sans Serif"
                context.fillText("circular DNA", centerX, centerY + 22)
            }

            function drawArc(context, centerX, centerY, radius, offset,
                             startAngle, endAngle, color) {
                context.strokeStyle = color
                context.lineWidth = Math.max(10, offset * 0.55)
                context.beginPath()
                context.arc(centerX, centerY, radius + offset, startAngle, endAngle)
                context.stroke()
            }

            onWidthChanged: requestPaint()
            onHeightChanged: requestPaint()
        }
    }

    Column {
        anchors.left: parent.left
        anchors.top: parent.top
        anchors.margins: 16
        spacing: 6

        Text {
            text: preview.selectedFeature
            color: preview.colors.ink
            font.pixelSize: 15
            font.weight: Font.DemiBold
        }
        Text {
            text: "Pinch or use controls to inspect"
            color: preview.colors.mutedInk
            font.pixelSize: 12
        }
    }

    Row {
        anchors.right: parent.right
        anchors.bottom: parent.bottom
        anchors.margins: 14
        spacing: 8

        Button {
            width: preview.colors.touchTarget
            height: preview.colors.touchTarget
            text: "−"
            Accessible.name: "Zoom out map"
            onClicked: preview.zoom = Math.max(0.7, preview.zoom - 0.15)
        }
        Button {
            width: preview.colors.touchTarget
            height: preview.colors.touchTarget
            text: "↺"
            Accessible.name: "Reset map zoom"
            onClicked: preview.zoom = 1.0
        }
        Button {
            width: preview.colors.touchTarget
            height: preview.colors.touchTarget
            text: "+"
            Accessible.name: "Zoom in map"
            onClicked: preview.zoom = Math.min(2.0, preview.zoom + 0.15)
        }
    }

    PinchHandler {
        target: plot
        minimumScale: 0.7
        maximumScale: 2.0
        onScaleChanged: preview.zoom = scale
    }

    TapHandler {
        onTapped: preview.featureActivated("Map selection at "
                                           + Math.round(point.position.x)
                                           + ", "
                                           + Math.round(point.position.y))
    }
}
