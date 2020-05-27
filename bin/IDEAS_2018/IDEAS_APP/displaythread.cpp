#include "displaythread.h"
#include <QFile>
#include <sstream>

void DisplayThread::start()
{
    int count = 0;
    QFile file("stdoutput_file");
    if (file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QTextStream in(&file);
        in.setCodec("UTF-8");
        QThread::sleep(2);
        for(int i =0; i<25 ; i++)
        {
            QThread::sleep(1);
            while(!in.atEnd())
            {
                //QThread::sleep(0.5);
                emit on_stdout(in.readLine());
                emit on_stdout_number(count*100/58);
                count++;
            }
        }
    }
}
