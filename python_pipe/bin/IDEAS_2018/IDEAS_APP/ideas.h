#ifndef IDEAS_H
#define IDEAS_H

#include <QMainWindow>
#include <QDebug>
#include "displaythread.h"

namespace Ui {
class IDEAS;
}

class IDEAS : public QMainWindow
{
    Q_OBJECT

public:
    explicit IDEAS(QWidget *parent = 0);
    ~IDEAS();


private slots:
    void on_goButton_clicked();

    //void convertStringToChar(std::string& str, char** array);

    void on_openFile_clicked();

    void init();

    void on_bedButton_clicked();

    void on_lineEdit_2_editingFinished();

    void on_lineEdit_3_editingFinished();

    void on_lineEdit_4_editingFinished();

    void on_lineEdit_5_editingFinished();

    void on_lineEdit_7_editingFinished();

    void on_lineEdit_8_editingFinished();

    void on_lineEdit_9_editingFinished();

    void on_lineEdit_10_editingFinished();

    void on_lineEdit_11_editingFinished();

    void changeColor();

    void createButtons(int, std::vector<std::string> labels, int number2, std::vector<std::string> colors);

    int getButtonNumber();

    void on_processChanges_clicked();

    void call_IDEAS(int argc, char* argv[]);

    void on_createTracksButton_clicked();

    void on_createHeatMapButton_clicked();

public slots:

    void stdout_print(QString name);

private:
    Ui::IDEAS *ui;
    DisplayThread displayThread;
};

#endif // IDEAS_H
