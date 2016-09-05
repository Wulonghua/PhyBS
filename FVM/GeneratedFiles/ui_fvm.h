/********************************************************************************
** Form generated from reading UI file 'fvm.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FVM_H
#define UI_FVM_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "RenderWidget.h"

QT_BEGIN_NAMESPACE

class Ui_FVMClass
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    RenderWidget *rWidget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *FVMClass)
    {
        if (FVMClass->objectName().isEmpty())
            FVMClass->setObjectName(QStringLiteral("FVMClass"));
        FVMClass->resize(743, 530);
        centralWidget = new QWidget(FVMClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        gridLayout->setContentsMargins(1, 1, 1, 1);
        rWidget = new RenderWidget(centralWidget);
        rWidget->setObjectName(QStringLiteral("rWidget"));

        gridLayout->addWidget(rWidget, 0, 0, 1, 1);

        FVMClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(FVMClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 743, 21));
        FVMClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(FVMClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        FVMClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(FVMClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        FVMClass->setStatusBar(statusBar);

        retranslateUi(FVMClass);

        QMetaObject::connectSlotsByName(FVMClass);
    } // setupUi

    void retranslateUi(QMainWindow *FVMClass)
    {
        FVMClass->setWindowTitle(QApplication::translate("FVMClass", "FVM", 0));
    } // retranslateUi

};

namespace Ui {
    class FVMClass: public Ui_FVMClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FVM_H
