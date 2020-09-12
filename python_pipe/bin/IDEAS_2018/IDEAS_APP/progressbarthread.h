#ifndef PROGRESSBARTHREAD_H
#define PROGRESSBARTHREAD_H


class progressBarThread
{
public:
    progressBarThread();

    void start();

signals:
    void progress(int);
};

#endif // PROGRESSBARTHREAD_H
