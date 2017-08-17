//    This is a GUI front-end for BayeScEnv
//    Copyright (C) 2014 Pierre de Villemereuil
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "gui.h"
#include "ui_gui.h"
#include <QDebug>
#include <QFileDialog>
#include <QProcess>
#include <QString>
#include <omp.h>

#if defined _WIN32 || defined _CYGWIN_
    #define windows 1
#else
    #define windows 0
#endif


GUI::GUI(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GUI)
{
    process = new QProcess(this);
    ui->setupUi(this);

    int threads;
    threads=omp_get_max_threads();
    ui->threadBox->setMaximum(threads);
    ui->threadBox->setMinimum(1);
    ui->threadBox->setValue(threads);
}

GUI::~GUI()
{
    process->kill();
    delete ui;
}

//-----------------------------------Misc functions

//Custom SLOT to connect STDOUT to the displayBox
void GUI::outlog()
{

    QString abc = process->readAllStandardOutput();
    emit outlogtext(abc);
    ui->displayBox->insertPlainText(abc);
}



//--------------------------------------Browsing

//Browse toward the input file
void GUI::on_browseInputButton_clicked()
{
    QString inputfile=QFileDialog::getOpenFileName(
                this,
                tr("Select the input file"),
                "",
                "All files (*.*)"
                );
    ui->inputEdit->setText(inputfile);
}

//Browse toward the discarded loci file
void GUI::on_browseDiscButton_clicked()
{
    QString discfile=QFileDialog::getOpenFileName(
                this,
                tr("Select the input file"),
                "",
                "All files (*.*)"
                );
    ui->discEdit->setText(discfile);
}

//Browse toward the environmental values file
void GUI::on_browseEnvButton_clicked()
{
    QString envfile=QFileDialog::getOpenFileName(
                this,
                tr("Select the input file"),
                "",
                "All files (*.*)"
                );
    ui->envEdit->setText(envfile);
}

//Browse toward the output folder
void GUI::on_browseOutfolderButton_clicked()
{
    QString outfolder = QFileDialog::getExistingDirectory(this,
                                                 tr("Select Output Directory"),
                                                 "",
                                                 QFileDialog::ShowDirsOnly
                                                 | QFileDialog::DontResolveSymlinks);
    ui->folderoutEdit->setText(outfolder);
}

//--------------------------------------Interactivity

void GUI::on_betaFisCheckBox_clicked()
{
    if(ui->betaFisCheckBox->isChecked()){
        ui->meanBetaBox->setEnabled(true);
        ui->sdBetaBox->setEnabled(true);
        ui->meanBetaLabel->setEnabled(true);
        ui->sdBetaLabel->setEnabled(true);
        ui->lowerFisBox->setEnabled(false);
        ui->upperFisBox->setEnabled(false);
        ui->lowerFisLabel->setEnabled(false);
        ui->upperFisLabel->setEnabled(false);
    } else {
        ui->meanBetaBox->setEnabled(false);
        ui->sdBetaBox->setEnabled(false);
        ui->meanBetaLabel->setEnabled(false);
        ui->sdBetaLabel->setEnabled(false);
        ui->lowerFisBox->setEnabled(true);
        ui->upperFisBox->setEnabled(true);
        ui->lowerFisLabel->setEnabled(true);
        ui->upperFisLabel->setEnabled(true);
    }
}

void GUI::on_FstatCheckBox_clicked()
{
   if(ui->FstatCheckBox->isChecked()){
        ui->envEdit->setEnabled(false);
        ui->envLabel->setEnabled(false);
        ui->UBox->setEnabled(false);
        ui->ULabel->setEnabled(false);
        ui->priorAlphaBox->setEnabled(false);
        ui->priorAlphaLabel->setEnabled(false);
        ui->prJumpBox->setEnabled(false);
        ui->prJumpLabel->setEnabled(false);
        ui->prPrefBox->setEnabled(false);
        ui->prPrefLabel->setEnabled(false);
    } else {
       ui->envEdit->setEnabled(true);
       ui->envLabel->setEnabled(true);
       ui->UBox->setEnabled(true);
       ui->ULabel->setEnabled(true);
       ui->priorAlphaBox->setEnabled(true);
       ui->priorAlphaLabel->setEnabled(true);
       ui->prJumpBox->setEnabled(true);
       ui->prJumpLabel->setEnabled(true);
       ui->prPrefBox->setEnabled(true);
       ui->prPrefLabel->setEnabled(true);
    }
}

//--------------------------------------Starting and quiting BayeScEnv

QString GUI::setup_call()
{
    QString call;
    if (windows)
        call="bayescenv.exe ";
    else
        call="./bayescenv ";

    //INPUT/OUTPUT section
    if (GUI::ui->inputEdit->text()!="")
    {
        call.append("\"");
        call.append(GUI::ui->inputEdit->text());
        call.append("\"");
    } else {
        return "ERROR";
    }
    if (GUI::ui->discEdit->text()!="") {
        call.append(" -d \"");
        call.append(GUI::ui->discEdit->text());
        call.append("\"");
    }
    if (GUI::ui->envEdit->text()!=""&&GUI::ui->envEdit->isEnabled()){
        call.append(" -env \"");
        call.append(GUI::ui->envEdit->text());
        call.append("\"");
    } else if (!ui->FstatCheckBox->isChecked()) {
        return "ERROR";
    }
    if (GUI::ui->folderoutEdit->text()!="") {
        call.append(" -od \"");
        call.append(GUI::ui->folderoutEdit->text());
        call.append("\"");
    }
    if (GUI::ui->nameoutEdit->text()!="") {
        call.append(" -o \"");
        call.append(GUI::ui->nameoutEdit->text());
        call.append("\"");
    }

    //OPTIONS section
    if (GUI::ui->SNPmatrixCheckBox->isChecked()) call.append(" -snp ");
    if (GUI::ui->FstatCheckBox->isChecked()) call.append(" -fstat ");
    if (GUI::ui->allTracesCheckBox->isChecked()) call.append(" -all_trace ");
    if (GUI::ui->outPilotCheckBox->isChecked()) call.append(" -out_pilot ");
    if (GUI::ui->outFreqCheckBox->isChecked()) call.append(" -out_freq ");


    //CHAIN section
    call.append(QString(" -threads %1").arg(ui->threadBox->value()));
    call.append(QString(" -burn %1").arg(ui->burninBox->value()));
    call.append(QString(" -thin %1").arg(ui->thinBox->value()));
    call.append(QString(" -n %1").arg(ui->sampleBox->value()));
    call.append(QString(" -nbp %1").arg(ui->nbpilotBox->value()));
    call.append(QString(" -pilot %1").arg(ui->lengthpilotBox->value()));

    //DOMINANT section
    if (GUI::ui->betaFisCheckBox->isChecked()) {
        call.append(" -beta_fis ");
        call.append(QString(" -m_fis %1").arg(ui->meanBetaBox->value()));
        call.append(QString(" -sd_fis %1").arg(ui->sdBetaBox->value()));
    } else {
        call.append(QString(" -lb_fis %1").arg(ui->lowerFisBox->value()));
        call.append(QString(" -hb_fis %1").arg(ui->upperFisBox->value()));
    }
    call.append(QString(" -aflp_pc %1").arg(ui->threshAFLPBox->value()));


    //GENERAL section
    call.append(QString(" -unif_pr %1").arg(ui->UBox->value()));
    call.append(QString(" -mean_alpha %1").arg(ui->priorAlphaBox->value()));
    call.append(QString(" -pr_jump %1").arg(ui->prJumpBox->value()));
    call.append(QString(" -pr_pref %1").arg(ui->prPrefBox->value()));

    return call;
}

//Start Button
void GUI::on_startButton_clicked()
{
    QString program=setup_call();
    if (program!="ERROR") {
        connect(process,SIGNAL(readyReadStandardOutput()),this,SLOT(outlog()));
        process->setProcessChannelMode(QProcess::MergedChannels);
        process->start(program);
        qDebug() << program;
        process->waitForReadyRead();
    } else {
        ui->displayBox->insertPlainText("ERROR: Some required fields (e.g. input or env file) are empty!\n");
    }
}

//Stop Button
void GUI::on_stopButton_clicked()
{
        process->kill();
        ui->displayBox->insertPlainText("\nStopped by the user!\n");
}
