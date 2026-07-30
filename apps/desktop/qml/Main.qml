import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import "theme" as Theme
import "components" as Components

ApplicationWindow {
    id: root

    visible: true
    width: 1240
    height: 800
    minimumWidth: 360
    minimumHeight: 620
    title: controller ? controller.appName : "OpenGene Studio"
    color: colors.canvas

    Theme.Theme { id: colors }

    property var controller: typeof appController === "undefined" ? null : appController
    property bool compact: width < 900
    property int currentWorkspace: 0
    property string currentWorkspaceName: workspaceNames[currentWorkspace]
    property var workspaceNames: ["Map", "Sequence", "Features", "Tools", "History"]
    property string previewStatus: "Native host ready — choose a workspace"

    function notify(message) {
        if (controller)
            controller.announce(message)
        else
            previewStatus = message
    }

    function chooseWorkspace(index) {
        currentWorkspace = index
        notify(workspaceNames[index] + " workspace selected")
    }

    header: ToolBar {
        id: topBar
        height: 64
        padding: 0

        background: Rectangle {
            color: colors.surface
            border.width: 1
            border.color: colors.border
        }

        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: root.compact ? 16 : 22
            anchors.rightMargin: root.compact ? 12 : 22
            spacing: 12

            Rectangle {
                Layout.preferredWidth: 34
                Layout.preferredHeight: 34
                radius: 11
                color: colors.accent

                Text {
                    anchors.centerIn: parent
                    text: "OG"
                    color: "white"
                    font.pixelSize: 12
                    font.weight: Font.Bold
                }
            }

            ColumnLayout {
                Layout.fillWidth: true
                spacing: 1

                Label {
                    text: "OpenGene Studio"
                    color: colors.ink
                    font.pixelSize: root.compact ? 16 : 18
                    font.weight: Font.DemiBold
                    elide: Text.ElideRight
                    Layout.fillWidth: true
                }
                Label {
                    visible: !root.compact
                    text: controller ? controller.coreStatus : "Native UI preview"
                    color: colors.mutedInk
                    font.pixelSize: 12
                    elide: Text.ElideRight
                    Layout.fillWidth: true
                }
            }

            ToolButton {
                implicitWidth: colors.touchTarget
                implicitHeight: colors.touchTarget
                text: "?"
                Accessible.name: "About native host"
                onClicked: aboutSheet.open()
            }
        }
    }

    footer: Components.TouchTabBar {
        visible: root.compact
        height: visible ? implicitHeight : 0
        width: parent ? parent.width : 0
        colors: colors
        currentIndex: root.currentWorkspace
        onDestinationRequested: root.chooseWorkspace(index)
    }

    RowLayout {
        anchors.fill: parent
        spacing: 0

        Components.NavigationRail {
            visible: !root.compact
            Layout.fillHeight: true
            colors: colors
            currentIndex: root.currentWorkspace
            onDestinationRequested: root.chooseWorkspace(index)
        }

        StackLayout {
            id: workspaceStack
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: root.currentWorkspace

            Flickable {
                clip: true
                contentWidth: width
                contentHeight: mapContent.implicitHeight + 40
                boundsBehavior: Flickable.DragOverBounds

                ColumnLayout {
                    id: mapContent
                    width: parent.width
                    spacing: 16
                    anchors.margins: root.compact ? 14 : colors.contentPadding

                    Label {
                        text: "Construct map"
                        color: colors.ink
                        font.pixelSize: 24
                        font.weight: Font.DemiBold
                    }
                    Label {
                        text: "Touch-first shell preview. Scientific map rendering will be bound to the portable core."
                        color: colors.mutedInk
                        font.pixelSize: 13
                        wrapMode: Text.WordWrap
                        Layout.fillWidth: true
                    }
                    Components.MapPreview {
                        Layout.fillWidth: true
                        Layout.preferredHeight: root.compact ? 380 : 500
                        colors: colors
                        onFeatureActivated: root.notify(label)
                    }
                    RowLayout {
                        Layout.fillWidth: true
                        spacing: 10
                        Button {
                            text: "Open project"
                            Layout.minimumHeight: colors.touchTarget
                            Layout.minimumWidth: 144
                            onClicked: root.notify("Project opening will be provided by the document controller")
                        }
                        Button {
                            text: "Import sequence"
                            Layout.minimumHeight: colors.touchTarget
                            Layout.minimumWidth: 154
                            onClicked: root.notify("Import workflow is a planned controller action")
                        }
                        Item { Layout.fillWidth: true }
                    }
                }
            }

            Flickable {
                clip: true
                contentWidth: width
                contentHeight: sequenceContent.implicitHeight + 40

                ColumnLayout {
                    id: sequenceContent
                    width: parent.width
                    spacing: 16
                    anchors.margins: root.compact ? 14 : colors.contentPadding

                    Label {
                        text: "Sequence"
                        color: colors.ink
                        font.pixelSize: 24
                        font.weight: Font.DemiBold
                    }
                    Components.WorkspaceCard {
                        Layout.fillWidth: true
                        colors: colors
                        title: "Base-level editor boundary"
                        detail: "A future document controller will supply visible base tiles, selection synchronization, undoable edit commands, and origin-aware coordinates."
                        actionLabel: "Preview command"
                        onActionRequested: root.notify("Sequence editing is intentionally not simulated in this host scaffold")
                    }
                    Rectangle {
                        Layout.fillWidth: true
                        Layout.preferredHeight: 190
                        radius: colors.cornerRadius
                        color: colors.surface
                        border.width: 1
                        border.color: colors.border
                        Text {
                            anchors.fill: parent
                            anchors.margins: 18
                            text: "001  ATGACGTTGACCTGATCGTACGATCGTAGCTAGCTAGCATCGATG\n051  CTAGCATCGATCGTACGATCGATGCTAGCTAGCATCGATCGATCG"
                            font.family: "Monospace"
                            font.pixelSize: root.compact ? 13 : 15
                            color: colors.ink
                            wrapMode: Text.WrapAnywhere
                        }
                    }
                }
            }

            Flickable {
                clip: true
                contentWidth: width
                contentHeight: featuresContent.implicitHeight + 40

                ColumnLayout {
                    id: featuresContent
                    width: parent.width
                    spacing: 12
                    anchors.margins: root.compact ? 14 : colors.contentPadding

                    Label {
                        text: "Features"
                        color: colors.ink
                        font.pixelSize: 24
                        font.weight: Font.DemiBold
                    }
                    Repeater {
                        model: [
                            ["Promoter", "1,022–1,184  •  forward strand"],
                            ["CDS", "1,225–2,789  •  forward strand"],
                            ["Origin of replication", "4,960–5,754  •  forward strand"]
                        ]
                        delegate: Components.WorkspaceCard {
                            required property var modelData
                            Layout.fillWidth: true
                            colors: colors
                            title: modelData[0]
                            detail: modelData[1] + "\nSample presentation only; annotations will come from the document model."
                            actionLabel: "Inspect"
                            onActionRequested: root.notify("Inspecting " + modelData[0])
                        }
                    }
                }
            }

            Flickable {
                clip: true
                contentWidth: width
                contentHeight: toolsContent.implicitHeight + 40

                ColumnLayout {
                    id: toolsContent
                    width: parent.width
                    spacing: 12
                    anchors.margins: root.compact ? 14 : colors.contentPadding

                    Label {
                        text: "Tools"
                        color: colors.ink
                        font.pixelSize: 24
                        font.weight: Font.DemiBold
                    }
                    Components.WorkspaceCard {
                        Layout.fillWidth: true
                        colors: colors
                        title: "Restriction analysis"
                        detail: "Connect the existing portable restriction package here; results must remain deterministic and document-scoped."
                        actionLabel: "Open tool"
                        onActionRequested: root.notify("Restriction analysis is awaiting controller wiring")
                    }
                    Components.WorkspaceCard {
                        Layout.fillWidth: true
                        colors: colors
                        title: "Primers and PCR"
                        detail: "Keep thermodynamic calculations in the portable engine and present queued jobs with explicit validation status."
                        actionLabel: "Open tool"
                        onActionRequested: root.notify("Primer and PCR tools are awaiting controller wiring")
                    }
                    Components.WorkspaceCard {
                        Layout.fillWidth: true
                        colors: colors
                        title: "Cloning operations"
                        detail: "Each future action returns a product and provenance record; this QML layer never mutates a sequence itself."
                        actionLabel: "Open tool"
                        onActionRequested: root.notify("Cloning operations are awaiting controller wiring")
                    }
                }
            }

            Flickable {
                clip: true
                contentWidth: width
                contentHeight: historyContent.implicitHeight + 40

                ColumnLayout {
                    id: historyContent
                    width: parent.width
                    spacing: 14
                    anchors.margins: root.compact ? 14 : colors.contentPadding

                    Label {
                        text: "Construct history"
                        color: colors.ink
                        font.pixelSize: 24
                        font.weight: Font.DemiBold
                    }
                    Components.WorkspaceCard {
                        Layout.fillWidth: true
                        colors: colors
                        title: "No operation selected"
                        detail: "The eventual history DAG will display immutable inputs, parameters, engine versions, warnings, hashes, and reproducibility status."
                        actionLabel: "Learn about provenance"
                        onActionRequested: root.notify("Provenance is represented by the portable history model")
                    }
                }
            }
        }
    }

    Rectangle {
        id: statusPill
        z: 10
        width: Math.min(root.width - 32, statusText.implicitWidth + 32)
        height: 38
        anchors.horizontalCenter: parent.horizontalCenter
        anchors.bottom: root.compact ? root.footer.top : parent.bottom
        anchors.bottomMargin: 14
        radius: 19
        color: colors.ink
        opacity: 0.94
        visible: true

        Text {
            id: statusText
            anchors.centerIn: parent
            width: parent.width - 24
            text: root.controller ? root.controller.statusMessage : root.previewStatus
            color: "white"
            font.pixelSize: 12
            elide: Text.ElideRight
            horizontalAlignment: Text.AlignHCenter
        }
    }

    Popup {
        id: aboutSheet
        modal: true
        focus: true
        anchors.centerIn: Overlay.overlay
        width: Math.min(root.width - 32, 420)
        padding: 20
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

        background: Rectangle {
            radius: colors.cornerRadius
            color: colors.surface
            border.width: 1
            border.color: colors.border
        }

        contentItem: ColumnLayout {
            spacing: 12
            Label {
                text: "Native host scaffold"
                color: colors.ink
                font.pixelSize: 20
                font.weight: Font.DemiBold
            }
            Label {
                text: "This is a clean-room Qt Quick shell. It does not simulate scientific results; a future controller will bind it to the portable OpenGene engines."
                color: colors.mutedInk
                font.pixelSize: 13
                wrapMode: Text.WordWrap
                Layout.fillWidth: true
            }
            Button {
                text: "Close"
                Layout.minimumHeight: colors.touchTarget
                Layout.minimumWidth: 112
                onClicked: aboutSheet.close()
            }
        }
    }
}
