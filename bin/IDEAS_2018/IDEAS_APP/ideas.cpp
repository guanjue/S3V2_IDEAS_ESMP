#include "ideas.h"
#include "ui_ideas.h"
#include "datastructure.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QColorDialog>
#include <QColor>
#include <QString>

#include <fstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>
#include <cstdlib>
#include <QProcess>
#include <qthread.h>
#include <QFuture>
#include <QtConcurrent/QtConcurrentRun>


#include "datastructure.h"
#include "genomicTensor.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


unsigned int rseed = 0;
bool hpass = false;
bool splitmerge = true;
bool SA = false;
bool gzip = true;
int bST = -1, bED = -1;
vector<string> imputelist;


IDEAS::IDEAS(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::IDEAS)
{
    ui->setupUi(this);
    //ui->progressBar->setValue(0);
    init();

    stdout = freopen("stdoutput_file", "w", stdout);
    //connect(ui->goButton, SIGNAL(readyReadStdOutput()), this, SLOT(updateText()));

}

IDEAS::~IDEAS()
{
    delete ui;
}

void IDEAS::init()
{
    ui->lcdNumber->setPalette(Qt::white);
    ui->seed_Slider->setMaximum(1000000);
    ui->lineEdit->setText("200");
    ui->lineEdit_2->setText("0");
    ui->lineEdit_3->setText("0");
    ui->lineEdit_4->setText("20");
    ui->lineEdit_5->setText("0");
    ui->lineEdit_6->setText("0");
    ui->lineEdit_7->setText("1");
    ui->lineEdit_8->setText("20");
    ui->lineEdit_9->setText("20");
    ui->lineEdit_10->setText("0.5");
    ui->lineEdit_11->setText("1000000");

}


/*................................. GO Button CLICK!.............................. */

void IDEAS::on_goButton_clicked()
{
    char * arg[22];
    arg[0] = "ideas";
    arg[1] = "t.input";
    arg[2] = "t.bed";
    arg[3] = "-log2";
    arg[4] = new char[10];
    strcpy(arg[4],(ui->lineEdit_2->text()).toStdString().c_str());   // log2
    arg[5] = "-G";
    arg[6] = new char[10];
    strcpy(arg[6],(ui->lineEdit_3->text()).toStdString().c_str());   // Maximum number of states to be inferred
    arg[7] = "-C";
    arg[8] = new char[10];
    strcpy(arg[8],(ui->lineEdit_4->text()).toStdString().c_str());   //Initial number of states
    arg[9] = "-P";
    arg[10] = new char[10];
    strcpy(arg[10],(ui->lineEdit_5->text()).toStdString().c_str());   //Initial number of states
    arg[11] = "-K";
    arg[12] = new char[10];
    strcpy(arg[12],(ui->lineEdit_6->text()).toStdString().c_str());   //Initial number of states
    arg[13] = "-A";
    arg[14] = new char[10];
    strcpy(arg[14],(ui->lineEdit_7->text()).toStdString().c_str());   //Initial number of states
    arg[15] = "-sample";
    arg[16] = new char[10];
    strcpy(arg[16],(ui->lineEdit_8->text()).toStdString().c_str());   //Initial number of states
    arg[17] = new char[10];
    strcpy(arg[17],(ui->lineEdit_9->text()).toStdString().c_str());   //Initial number of states
    arg[18] = "-minerr";
    arg[19] = new char[10];
    strcpy(arg[19],(ui->lineEdit_10->text()).toStdString().c_str());   //Initial number of states
    arg[20] = "-maxerr";
    arg[21] = new char[10];
    strcpy(arg[21],(ui->lineEdit_11->text()).toStdString().c_str());   //Initial number of states

   // call_IDEAS(22,(char**)arg);
   // IDEAS ideas;

    QFuture<void> t2 = QtConcurrent::run(this,&IDEAS::call_IDEAS,22,(char**)arg);

    //connect(&displayThread, SIGNAL(on_stdout_number(int)),ui->progressBar,SLOT(setValue(int)));
    connect(&displayThread, &DisplayThread::on_stdout,this,&IDEAS::stdout_print);

    //delete [] arg;

    //int numberOfButtons = getButtonNumber();

    QFuture<void> t1 = QtConcurrent::run(&this->displayThread, &DisplayThread::start);
    //QFuture<void> t3 = QtConcurrent::run(&this->prog, &progressBarThread::start);
    //cout << "t2 finished" << endl;
    //t1.waitForFinished();

    //t2.waitForFinished();
}
/*{
    if(t2.isFinished())
        {

        //.......................Create Buttons for buttons below the picture..............

            std::ifstream file("t.input");
            string firstWord;
            std::vector<std::string> words;

            string a,b,c;

            //string word;
            int count = 1;
            file >> a >> b >> c;
            words.push_back(b);
            firstWord = b;

            while(!file.eof())
            {
                //file >> word;

                //if (word == firstWord)
                file >> a >> b >> c;
                if (b == firstWord)
                    break;
                words.push_back(b);
                count++;
            }

            file.close();

            system("/usr/bin/Rscript /home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/run.R"); // Call R Script to create the heatmap

            QPixmap pix("/home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/Rplot.png");
            ui->ImageLabel->setPixmap(pix);

            //.......................Create Buttons for buttons beside the picture................

            std::ifstream statecolorfile("statecolor.txt");
            std::vector<std::string> colors;

            string a1,b1;
            int count1=0;

            while(!statecolorfile.eof())
            {
                statecolorfile >> a1 >> b1;
                colors.push_back(a1);
                count1++;
            }

            statecolorfile.close();

            createButtons(count, words, count1-1, colors);    // Call createButtons function to create both set of buttons

        }
    }
*/
/*................................. Open Button CLICK!.............................. */

void IDEAS::on_openFile_clicked()
{
    QString filename = QFileDialog::getOpenFileName();

    QFile file(filename);
     if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
         return;

    QString content = file.readAll();

    file.close();
}

/*................................. Load Bed file CLICK!.............................. */

void IDEAS::on_bedButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName();

    QFile file(filename);
     if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
         return;

    QString content = file.readAll();

    file.close();
}


/*................................. Set Default Values................................. */

void IDEAS::on_lineEdit_2_editingFinished()
{
    QString text = ui->lineEdit_2->text();
    bool ok;
    double value = text.toDouble(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Use log2(x+number) transformation"),
                tr("Value is not a number") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Use log2(x+number) transformation"),
                tr("Value should be greater than or equal to zero") );
    }

}

void IDEAS::on_lineEdit_3_editingFinished()
{
    QString text = ui->lineEdit_3->text();
    bool ok;
    int value = text.toInt(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Maximum number of states to be inferred"),
                tr("Value is not an integer") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Maximum number of states to be inferred"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_4_editingFinished()
{
    QString text = ui->lineEdit_4->text();
    bool ok;
    int value = text.toInt(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Initial number of states"),
                tr("Value is not an integer") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Initial number of states"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_5_editingFinished()
{
    QString text = ui->lineEdit_5->text();
    bool ok;
    double value = text.toDouble(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Maximum number of position classes to be inferred"),
                tr("Value is not a number") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Maximum number of position classes to be inferred"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_7_editingFinished()
{
    QString text = ui->lineEdit_7->text();
    bool ok;
    double value = text.toDouble(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Prior concentration"),
                tr("Value is not a number") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Prior concentration"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_8_editingFinished()
{
    QString text = ui->lineEdit_8->text();
    bool ok;
    int value = text.toInt(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Number of burnin steps"),
                tr("Value is not an integer") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Number of burnin steps"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_9_editingFinished()
{
    QString text = ui->lineEdit_9->text();
    bool ok;
    int value = text.toInt(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Number of maximization steps"),
                tr("Value is not an integer") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Number of maximization steps"),
                tr("Value should be greater than or equal to zero") );
    }
}

void IDEAS::on_lineEdit_10_editingFinished()
{
    QString text = ui->lineEdit_10->text();
    bool ok;
    double value = text.toDouble(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Minimum standard deviation for the emission Gaussian distribution"),
                tr("Value is not a number") );
    }

    if (value <= 0)
    {
        QMessageBox::information(
                this,
                tr("Minimum standard deviation for the emission Gaussian distribution"),
                tr("Value should be greater than zero") );
    }
}

void IDEAS::on_lineEdit_11_editingFinished()
{
    QString text = ui->lineEdit_11->text();
    bool ok;
    double value = text.toDouble(&ok);
    if (!ok)
    {
        QMessageBox::information(
                this,
                tr("Maximum standard deviation for the emission Gaussian distribution"),
                tr("Value is not a number") );
    }

    if (value < 0)
    {
        QMessageBox::information(
                this,
                tr("Maximum standard deviation for the emission Gaussian distribution"),
                tr("Value should be greater than zero") );
    }
}

/*................................. Change button color................................. */

void IDEAS::changeColor()
{
    QColor j = ((QPushButton*)sender())->backgroundRole();

    if(j.red() == 0 && j.green()==0 && j.blue()==1)
    {
        std::ofstream outfile1;
        std::ofstream outfile2;
        outfile1.open("mc.txt",std::ofstream::out | std::ofstream::app);
        outfile2.open("sc.txt",std::ofstream::out | std::ofstream::app);
        QColor color = QColorDialog::getColor(Qt::white, this);
        if(color.isValid())
        {
            QString hexColor = color.name();
            //printf("");
            //ui->button[i]->setPalette(color);
            QString red = QString::number(color.red());
            QString green = QString::number(color.green());
            QString blue = QString::number(color.blue());
            QString rgb = "("+red+","+green+","+blue+")";
            ((QPushButton*)sender())->setStyleSheet("background-color:rgb"+rgb);
            QString obj = ((QPushButton*)sender())->objectName();
            if(obj.toInt()<50)
            {
                outfile1 << obj.toStdString() +" "+ red.toStdString()+" "+green.toStdString()+" "+blue.toStdString()+"\n";
                outfile1.close();
            }
            else
            {
                outfile2 << obj.toStdString() +" "+ red.toStdString()+","+green.toStdString()+","+blue.toStdString()+" "+hexColor.toStdString()+"\n";
                outfile2.close();


                system("/usr/bin/Rscript /home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/run.R NULL sc.txt"); // Call R Script to create the heatmap

                QPixmap pix("/home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/Rplot.png");
                ui->ImageLabel->setPixmap(pix);

                //std::ofstream ofs;
                //ofs.open("sc.txt", std::ofstream::out | std::ofstream::trunc);
                //ofs.close();

            }

            //colorVec[obj.toInt()] = rgb.toStdString();
        }
    }

    else
        ((QPushButton*)sender())->setPalette(Qt::blue);
}

/*................................. Create Buttons ...................................... */

void IDEAS::createButtons(int number, std::vector<std::string> labels, int number2, std::vector<std::string> colors)
{
    QPushButton *button[number];
    QLabel *label[number];
    for(int i = 0; i< number; i++)
    {
      button[i] = new QPushButton;
      button[i]->setMaximumWidth(15);
      button[i]->setMinimumWidth(10);
      button[i]->setMaximumHeight(100);
      button[i]->setMinimumHeight(50);
      button[i]->setObjectName(QString::number(i+1));
      // configure your button with the common settings here
      ui->horizontalLayout->addWidget(button[i]);
      //button[i]->setPalette(Qt::blue);

      QString buttonColor = QString::fromStdString(colors[(i)]);
      button[i]->setStyleSheet("background-color:rgb("+buttonColor+")"); ////// Uncomment for coloring the state color buttons

      connect(button[i],SIGNAL(clicked()),this,SLOT(changeColor()));
      //connect(button[i],SIGNAL(),this, SLOT(changeButtonColor()));
    }

    QPushButton *buttonHorizontal[number2];
    //ui->verticalLayout_2->set
    for(int i = 0; i< number2; i++)
    {
      buttonHorizontal[i] = new QPushButton;
      buttonHorizontal[i]->setMaximumWidth(50);
      buttonHorizontal[i]->setMinimumWidth(15);
      buttonHorizontal[i]->setMaximumHeight(10);
      buttonHorizontal[i]->setMinimumHeight(5);
      buttonHorizontal[i]->setObjectName(QString::number(number2-i+50));
      ui->verticalLayout_2->addWidget(buttonHorizontal[i]);
      //if(number2 > 10)
      ui->verticalLayout_2->setSpacing((370/number2)-15);
      //else
      //    ui->verticalLayout_2->setSpacing(number2*1.5);
      //buttonHorizontal[i]->setPalette(Qt::blue);
      //QString text = QString::fromStdString(colors[(number2-1)-i]);
      //buttonHorizontal[i]->setStyleSheet("background-color:rgb("+text+")"); ////// Uncomment for coloring the state color buttons

      connect(buttonHorizontal[i],SIGNAL(clicked()),this,SLOT(changeColor()));
    }

    /*
    for(int j =0; j< number ; j++)
    {
        label[j] = new QLabel;
        ui->horizontalLayout_2->addWidget(label[j]);
        QString text = QString::fromStdString(labels[j]);
        label[j]->setText(text);
    }
    */
}


/*................................. Get Number of Buttons to be created.............................. */
int IDEAS::getButtonNumber()
{
    std::ifstream file("t.input");
    string firstWord;

    string a,b,c;
    file >> a >> b >> c;
    string word;
    int count = 0;
    file >> firstWord;

    while(!file.eof())
    {
        file >> word;
        if (word == firstWord)
            count++;
    }

    file.close();
    return count;
}

/*
void IDEAS::on_refresh_mc_clicked()
{
    std::ofstream ofs;
    ofs.open("mc.txt", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
}

void IDEAS::on_refresh_sc_clicked()
{
    std::remove("sc.txt");

}
*/

void IDEAS::on_processChanges_clicked()
{
    system("/usr/bin/Rscript /home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/run.R mc.txt NULL"); // Call R Script to create the heatmap

    QPixmap pix("/home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/Rplot.png");
    ui->ImageLabel->setPixmap(pix);

    std::ofstream ofs;
    ofs.open("mc.txt", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
}


void IDEAS::stdout_print(QString name)
{
    ui->outputTextEdit->append(name);
}

void IDEAS::call_IDEAS(int argc, char* argv[])
{    
    int i, j;
    int minCut = 5, admixN = 30, burnin = 50, mcmc = 50, maxHapK = 0, maxGG = 0, maxPos = 0, thread = 1;
    bool sqc = false, lik = false, samplemaximum = true, outputproportion = false;
    double log2 = -1;
    bool add2 = false, ind = false, indind = false, nb = false, norm = false;
    double error = 0.01, A = 0, recr = 100., heteroVh = 1., minerr = 0.5, maxerr = 100000000.;
    char const *output = argv[1], *fixPop = NULL, *input = NULL, *fbed = NULL, *fcov = NULL, *fcovbed = NULL, *fparam = NULL, *fmixpara = NULL, *parafile = NULL, *statefile = NULL, *clusterfile = NULL, *para0file = NULL, *state0file = NULL, *cluster0file = NULL, *profile0file = NULL;
    int startK=0, fixC=0;
    vector<string> fS, fA;
    vector<vector<int> > remove;

    i = 1;
    input = argv[1];
    if(argc > 2 && argv[2][0] != '-')
    {	fbed = argv[2];
        i = 2;
    }

    for(i = i + 1; i < argc; i++)
    {	if(strcmp(argv[i], "-rseed") == 0)
        {	rseed = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-sa") == 0)
        {	SA = true;
        }
        else if(strcmp(argv[i], "-o") == 0)
        {	output = argv[i + 1];
            i++;
        }
        else if(strcmp(argv[i], "-P") == 0)
        {	maxPos = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-K") == 0)
        {	maxHapK = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-G") == 0)
        {	maxGG = atoi(argv[i + 1]);
            startK = max(startK, maxGG);
            i++;
        }
        else if(strcmp(argv[i], "-c") == 0)
        {	minCut = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-err") == 0)
        {	error = atof(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-inv") == 0)
        {	bST = atoi(argv[i + 1]);
            bED = atoi(argv[i + 2]);
            i += 2;
        }
        else if(strcmp(argv[i], "-nogz") == 0)
        {	gzip = false;
        }
        else if(strcmp(argv[i], "-impute") == 0)
        {	int l = (int)strlen(argv[i + 1]);
            int ii, jj;
            jj = 0;
            for(ii = 0; ii < l; ii++)
            {	if(argv[i + 1][ii] == ',')
                {	argv[i + 1][ii] = 0;
                    imputelist.push_back(&argv[i + 1][jj]);
                    jj = ii + 1;
                    argv[i + 1][ii] = ',';
                }
            }
            if(jj < ii) imputelist.push_back(&argv[i + 1][jj]);
            //for(ii = 0; ii < (int)imputelist.size(); ii++)
            //	printf("%s\n", imputelist[ii].c_str());
            i++;
        }
        else if(strcmp(argv[i], "-fixPop") == 0)
        {	fixPop = argv[i + 1];
            i++;
        }
        else if(strcmp(argv[i], "-fixHap") == 0)
        {	if((int)fS.size() != (int)fA.size()) printf("-fixHap option cannot be used together with -fixAllele option.\n"),exit(0);
            fS.push_back(argv[i + 1]);
            fA.push_back(argv[i + 2]);
            i += 2;
        }
        else if(strcmp(argv[i], "-fixAllele") == 0)
        {	if((int)fS.size() > 0) printf("-fixAllele option cannot be used together with -fixHap option.\n"),exit(0);
            fA.push_back(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-fixC") == 0)
        {	fixC = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-admixn") == 0)
        {	admixN = atoi(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-sqc") == 0)
            sqc = true;
        else if(strcmp(argv[i], "-sample") == 0)
        {	mcmc = atoi(argv[i + 2]);
            burnin = atoi(argv[i + 1]);
            i += 2;
        }
        else if(strcmp(argv[i], "-lik") == 0)
            lik = true;
        else if(strcmp(argv[i], "-A") == 0)
        {	A = atof(argv[i + 1]);
            i ++;
        }
        else if(strcmp(argv[i], "-z") == 0)
        {	fcov = argv[i + 1];
            i++;
            if(i + 1 < argc && argv[i + 1][0] != '-')
            {	fcovbed = argv[i + 1];
                i++;
            }
        }
        else if(strcmp(argv[i], "-max") == 0)
        {	samplemaximum = true;
        }
        else if(strcmp(argv[i], "-outputprop") == 0)
        {	outputproportion = true;
        }
        else if(strcmp(argv[i], "-rec") == 0)
        {	recr = atof(argv[i + 1]) * 1000.;
            i++;
        }
        else if(strcmp(argv[i], "-log2") == 0)
        {	log2 = 1.;
            if(argc > i + 1 && argv[i + 1][0] >= 48 && argv[i + 1][0] < 58)
            {	log2 = atof(argv[i + 1]);
                i++;
                //printf("log2 constant = %f\n", log2);
            }
        }
        else if(strcmp(argv[i], "-minerr") == 0)
        {	minerr = atof(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-maxerr") == 0)
        {	maxerr = atof(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-add2") == 0) add2 = true;
        else if(strcmp(argv[i], "-param") == 0)
        {	fparam = argv[i + 1];
            i++;
        }
        else if(strcmp(argv[i], "-mixpara") == 0)
        {	fmixpara = argv[i + 1];
            i++;
        }
        else if(strcmp(argv[i], "-startpara") == 0) //start with specified gauss parameters
        {	parafile = argv[i + 1];
            i++;
        }
        else if(strcmp(argv[i], "-prevrun") == 0) //restore segments from previous run on the same data and then continue to run
        {	statefile = argv[i + 1];
            clusterfile = argv[i + 2];
            i += 2;
        }
        else if(strcmp(argv[i], "-otherpara") == 0) //use other gauss parameters as priors
        {	para0file = argv[i + 1];
            i++;
            if(argc > i + 1 && argv[i + 1][0] != '-')
            {	profile0file = argv[i + 1];
                i++;
            }
        }
        else if(strcmp(argv[i], "-otherstate") == 0) //use other segments as priors
        {	state0file = argv[i + 1];
            cluster0file = argv[i + 2];
            i += 2;
        }
        else if(strcmp(argv[i], "-norm") == 0)
        {	norm = true;
        }
        else if(strcmp(argv[i], "-C") == 0)
        {	startK=atoi(argv[i+1]);
            i++;
        }
        else if(strcmp(argv[i], "-ind") == 0)
        {	ind = true;
        }
        else if(strcmp(argv[i], "-indind") == 0)
        {	indind = true;
        }
        else if(strcmp(argv[i], "-h") == 0)
        {	heteroVh = atof(argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-hp") == 0)
        {	hpass = true;
        }
        else if(strcmp(argv[i], "-nosm") == 0)
        {	splitmerge = false;
        }
        else if(strcmp(argv[i], "-remove") == 0)
        {	i++;
            int ind = atoi(argv[i]), ll = (int)strlen(argv[i]) - 1;
            vector<int> row;
            for(j = 1; j < ll; j++) if(argv[i][j] == ':') break;
            do {
                row.push_back(atoi(&argv[i][j + 1]));
                for(j = j + 1; j < ll; j++) if(argv[i][j] == ',') break;
            } while(j < ll);
            if((int)remove.size() <= ind) remove.resize(ind + 1, vector<int>());
            remove[ind] = row;
        }
        else if(strcmp(argv[i], "-nb") == 0)
        {	nb = true;
        }
        else if(strcmp(argv[i], "-thread") == 0)
        {	thread = atoi(argv[i + 1]);
            i++;
        }
    }

    time_t tst, ted;
    time(&tst);
    genomicTensor gtm(rseed);
    gtm.outputproportion = outputproportion;
    gtm.log2 = log2;
    gtm.norm = norm;
    if(fixC > 0) startK = fixC;
    gtm.run(input, fbed, fcov, fcovbed, burnin, mcmc, thread, maxHapK, maxGG, maxPos, A, recr, heteroVh, samplemaximum, output, sqc, fparam, add2, minerr, maxerr, ind, indind, startK, fixC, remove, nb, fmixpara, parafile, statefile, clusterfile, para0file, state0file, cluster0file, profile0file);

    time(&ted);
    //printf("total time = %dsec\n", (int)(ted - tst));

    //return 0
}


void IDEAS::on_createTracksButton_clicked()
{
    stdout = freopen("stdoutput_file", "w", stdout);

    system("/usr/bin/Rscript /home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/createTracks.R"); // Call R Script to create the heatmap

    QFuture<void> t1 = QtConcurrent::run(&this->displayThread, &DisplayThread::start);

}

void IDEAS::on_createHeatMapButton_clicked()
{
    //.......................Create Buttons for buttons below the picture..............

        std::ifstream file("t.input");
        string firstWord;
        std::vector<std::string> words;

        string a,b,c;

        //string word;
        int count = 1;
        file >> a >> b >> c;
        words.push_back(b);
        firstWord = b;

        while(!file.eof())
        {
            //file >> word;

            //if (word == firstWord)
            file >> a >> b >> c;
            if (b == firstWord)
                break;
            words.push_back(b);
            count++;
        }

        file.close();

        system("/usr/bin/Rscript /home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/run.R"); // Call R Script to create the heatmap

        QPixmap pix("/home/sachin/Downloads/IDEAS_Test/build-IDEAS_APP-Desktop_Qt_5_7_0_GCC_64bit-Debug/Rplot.png");
        ui->ImageLabel->setPixmap(pix);

        //.......................Create Buttons for buttons beside the picture................

        std::ifstream markcolorfile("markcolor.txt");
        std::vector<std::string> colors;

        string a1,b1;
        int count1=0;

        while(!markcolorfile.eof())
        {
            markcolorfile >> a1;
            colors.push_back(a1);
            count1++;
        }

        markcolorfile.close();

        createButtons(count, words, count1-1, colors);    // Call createButtons function to create both set of buttons

        fclose(stdout);

}
