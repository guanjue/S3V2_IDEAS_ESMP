#ifndef DISPLAYTHREAD_H
#define DISPLAYTHREAD_H

#include <QObject>
#include <QDebug>
#include <QThread>
#include <QFile>
#include <iostream>

class DisplayThread : public QObject
{
    Q_OBJECT
public:

    void start();

signals:
    void on_stdout(QString);
    void on_stdout_number(int);

 };


#endif // DISPLAYTHREAD_H
