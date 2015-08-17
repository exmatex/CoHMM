#include "RedisWrapper.hpp"

#include <cstdarg>
#include <iostream>

#include "hiredis.h"

RedisWrapper::RedisWrapper()
{
    this->hostSet = false;
}

RedisWrapper::~RedisWrapper()
{
    redisFree(this->rContext);
}

RedisWrapper& RedisWrapper::getContext(const char * redisHost)
{
    static RedisWrapper context;
    ///TODO: is it safe to check this variable outside a lock?
    //context.rMutex.lock();
    if(context.hostSet == false)
    {
        context.setHost(redisHost);
    }
    //context.rMutex.unlock();

    return context;
}

void RedisWrapper::setHost(const char * redisHost)
{
    rMutex.lock();
    //Possible that task spun before other thread set host
    if(this->hostSet == false)
    {
        this->rContext = redisConnect(redisHost, 6379);
        if(this->rContext != nullptr && this->rContext->err)
        {
            std::cerr << "Error: " << this->rContext->errstr << std::endl;
        }
        this->hostSet = true;
    }
    rMutex.unlock();
}


redisReply * RedisWrapper::redisCommand(const char *format, ...)
{
    va_list ap;
    redisReply * reply = nullptr;
    va_start(ap, format);
    rMutex.lock();
    reply = (redisReply *) redisvCommand(this->rContext, format, ap);
    rMutex.unlock();
    va_end(ap);
    return reply;
}
