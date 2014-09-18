#ifndef REDISBUCKETS_HPP
#define REDISBUCKETS_HPP

#include <hiredis.h>
#include <set>
#include <map>
#include <string>
#include <vector>

char * buildKey(double w0[7], char * tag, int keyDigits);
void putData(double w0[7], double f0[7], double g0[7], char * tag, redisContext *redis, int keyDigits);
void getBucket(char * key, redisContext *redis, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut);
void getBucket(double w0[7], char * tag, redisContext *redis, int keyDigits, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut);
void getSortedSubBucket(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut);
std::set<char *> * rebuildDB(std::set<char *> * keyDict, char * tag, redisContext *redis, int newKeyDigits);
void getSortedSubBucketNearZero(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut, double zeroThresh);
void getCachedSortedSubBucketNearZero(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut, double zeroThresh, std::map<std::string, std::vector<char *> > *dbCache);
void writeMessage(char * message, redisContext *redis);
#endif

