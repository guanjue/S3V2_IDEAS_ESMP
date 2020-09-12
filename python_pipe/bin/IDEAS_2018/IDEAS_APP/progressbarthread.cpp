#include "progressbarthread.h"

progressBarThread::progressBarThread()
{

}

void progressBarThread::start()
{
    int count = 0;

    std::ifstream file("stdoutput_file");

    string word = "itern";
    while(file >> word)
    {
        count++;
        emit progress((count*100) / 40);

    }

    file.close();
}
