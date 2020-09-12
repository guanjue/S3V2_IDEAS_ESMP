#include "ideas.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    IDEAS w;
    w.show();

    return a.exec();
}
