/** redis database utilities
 * calc buckets
 * sort
 * (un)pack
 * **/

#include "redisBuckets.hpp"

#include <hiredis.h>
#include <vector>
#include <set>
#include <cstring>
#include <cassert>
#include <list>
#include <cstdio>
#include <cmath>
#include "2DKriging.hpp"

struct SortableRedisEntry_s
{
	double w[7];
	double f[7];
	double g[7];
	double distance;
};

bool SortableRedisEntry_c(SortableRedisEntry_s &lhs, SortableRedisEntry_s &rhs)
{
	return lhs.distance < rhs.distance;
}

void redisWrite_fields(Node* nodes, Input in, redisContext* redis, int cur_step)
{
	redisReply * reply;
    int grid_size = in.dim_x*in.dim_y;
	
	//Buld the keys
	char nKey[1024];
	char stepKey[1024];
	sprintf(nKey, "macro:%d:nodes:", cur_step);
	sprintf(stepKey, "macro:step:");
	
	//Write W, F, and G
	reply = (redisReply *)redisCommand(redis, "HMSET %s stat %b", nKey, nodes, sizeof(Node)*grid_size);
    //printf("Save nodes %s\n", reply->str);//that shows +OK
	freeReplyObject(reply);

	//Update the current timestep in the database.
	reply = (redisReply *)redisCommand(redis, "HMSet %s stat %b", stepKey, &cur_step, sizeof(int));
    //put input with curr grad_threshold in db
	freeReplyObject(reply);
    Save_Input tmp = in;
	reply = (redisReply *)redisCommand(redis, "HMSet %s stat1 %b", stepKey, &tmp, sizeof(Save_Input));
    //printf("Save step %s\n", reply->str);//that shows +OK
	freeReplyObject(reply);

}

int redisRead_fields(Node* nodes, Input* in, redisContext * redis)
{
    int saved_step = 0;
	redisReply * reply;
	char stepKey[1024];
	char nKey[1024];
    //get last saved integration step
	sprintf(stepKey, "macro:step:");
	reply = (redisReply *)redisCommand(redis, "HMGet %s stat", stepKey);
    int* save_buffer = (int*)reply->element[0]->str;
    if(save_buffer){
	    memcpy(&saved_step, &save_buffer[0], sizeof(int));
        //printf("savedstep: %d\n", saved_step);
	    sprintf(nKey, "macro:%d:nodes:", saved_step);
	    freeReplyObject(reply);
	    reply = (redisReply *)redisCommand(redis, "HMGet %s stat1", stepKey);
        Save_Input* in_buffer = (Save_Input*)reply->element[0]->str;
        //Save_Input tmp;
	    //memcpy(&tmp, &in_buffer[0], sizeof(Save_Input));
	    memcpy(in, &in_buffer[0], sizeof(Save_Input));
        int grid_size = in->dim_x*in->dim_y;
        //printf("gradthresh %lf \n", in->grad_threshold);
	    freeReplyObject(reply);
        //load node field
	    reply = (redisReply *)redisCommand(redis, "HMGet %s stat", nKey);
        Node* node_buffer = (Node*)reply->element[0]->str;
	    memcpy(nodes, &node_buffer[0], sizeof(Node)*grid_size);
        //printf("nodes %lf\n", nodes[1].f.f[6]);
	    freeReplyObject(reply);

        //CkExit();
        printf("Loaded field of integration step %d\n", saved_step);
	return saved_step;
    }
    else{
        printf("No field loaded\n");
	    return -1;
    }
}

void redisDel_fields(redisContext* redis, int cur_step){

  redisReply * reply;
  if(cur_step>5){
	char nKey[1024];
	sprintf(nKey, "macro:%d:nodes:", cur_step-5);
	//Delete field
	reply = (redisReply *)redisCommand(redis, "HDEL %s stat", nKey);
    //printf("Del fields %s\n", reply->str);//that shows +OK
	freeReplyObject(reply);

  }
}

void buildKey(char* key, double w0[7], char * tag, int keyDigits)
{
	char fBuff[1024];
	//Build format
	sprintf(fBuff, "%%.%dlf:%%.%dlf:%%.%dlf:%%.%dlf:%%.%dlf:%%.%dlf:%%.%dlf:%%s:", keyDigits, keyDigits, keyDigits, keyDigits, keyDigits, keyDigits, keyDigits);
	//Use format
	sprintf(key, fBuff, w0[0], w0[1], w0[2], w0[3], w0[4], w0[5], w0[6], tag);
	//Return the key
}

void packValue(double w[7], double f[7], double g[7], double* retBuffer)
{
	//double * retBuffer = new double[21];
	//Memcpy all that stuff
	memcpy(&retBuffer[0], w, sizeof(double)*7);
	memcpy(&retBuffer[7], f, sizeof(double)*7);
	memcpy(&retBuffer[14], g, sizeof(double)*7);

	//Return it, switching to chars/bytes because that be what I likes
	//return (char *) retBuffer;
}

void unpackValue(char * value, double w[7], double f[7], double g[7])
{
	double * srcBuffer = (double *)value;
	//Memcpy all that stuff
	memcpy(w, &srcBuffer[0], sizeof(double)*7);
	memcpy(f, &srcBuffer[7], sizeof(double)*7);
	memcpy(g, &srcBuffer[14], sizeof(double)*7);

	return;
}

void putData(double w0[7], double f0[7], double g0[7], char * tag, redisContext *redis, int keyDigits)
{
	//Build Key
	char * key = new char[1024];
	buildKey(key, w0, tag, keyDigits);
	//Build Value
	double * retBuffer = new double[21];
	packValue(w0, f0, g0, retBuffer);
	//Redis the key/value pair
	redisReply *reply;
	reply = (redisReply *) redisCommand(redis, "SADD %s %b", key, retBuffer, sizeof(double)*21);
	freeReplyObject(reply);
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO SAVED:DATA:w:%e:%e:%e:%e:%e:%e:%e", w0[0], w0[1], w0[2], w0[3], w0[4], w0[5], w0[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO SAVED:DATA:f:%e:%e:%e:%e:%e:%e:%e", f0[0], f0[1], f0[2], f0[3], f0[4], f0[5], f0[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO SAVED:DATA:g:%e:%e:%e:%e:%e:%e:%e", g0[0], g0[1], g0[2], g0[3], g0[4], g0[5], g0[6]);
        freeReplyObject(traceRep);
#endif
  //delete value;
  delete[] key;
  delete[] retBuffer;
	return;
}

void getBucket(char * key, redisContext *redis, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut)
{
	//Get values	
	redisReply *reply;
	reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", key);
	assert(reply->type == REDIS_REPLY_ARRAY);
	//Write values to vectors
	int nVals = reply->elements;
	if(nVals == 0)
	{
		freeReplyObject(reply);
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO MISS:%s", key);
        freeReplyObject(traceRep);
#endif
		return;
	}
	else
	{
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO HIT:%s:values:%i", key,nVals);
        freeReplyObject(traceRep);
#endif
		for(int i=0; i<nVals; ++i)
		{
			double * writeW = new double[7];
			double * writeF = new double[7];
			double * writeG = new double[7];

			char *val = (reply->element)[i]->str;
			unpackValue(val, writeW, writeF, writeG);
#ifdef TRACE
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:w:%e:%e:%e:%e:%e:%e:%e", i, writeW[0], writeW[1], writeW[2], writeW[3], writeW[4], writeW[5], writeW[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:f:%e:%e:%e:%e:%e:%e:%e", i, writeF[0], writeF[1], writeF[2], writeF[3], writeF[4], writeF[5], writeF[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:g:%e:%e:%e:%e:%e:%e:%e", i, writeG[0], writeG[1], writeG[2], writeG[3], writeG[4], writeG[5], writeG[6]);
        freeReplyObject(traceRep);
#endif
			wOut->push_back(writeW);
			fOut->push_back(writeF);
			gOut->push_back(writeG);
            delete[] writeW;
            delete[] writeF;
            delete[] writeG;
		}
		freeReplyObject(reply);
		return;
	}
}

void getBucket(double w0[7], char * tag, redisContext *redis, int keyDigits, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut)
{
	//Build Key
	char * key = new char[1024];
	buildKey(key, w0, tag, keyDigits);
	//Use the key
	getBucket(key, redis, wOut, fOut, gOut);
  delete[] key;
	return;	
}

void getSortedSubBucket(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut)
{
	//Build Key
	char * key = new char[1024];
	buildKey(key, w0, tag, keyDigits);
	//Get values
	redisReply *reply;
	reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", key);
	if(reply->type != REDIS_REPLY_ARRAY)
	{
		freeReplyObject(reply);
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO MISS:%s", key);
        freeReplyObject(traceRep);
#endif
		return;
	}
	//Write values to set to sort
	std::list<SortableRedisEntry_s> sortSet;
	int nVals = reply->elements;
	if(nVals == 0)
	{
		freeReplyObject(reply);
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO MISS:%s", key);
        freeReplyObject(traceRep);
#endif
		return;
	}
	else
	{
#ifdef TRACE
        redisReply * traceRep;
        traceRep = (redisReply *)redisCommand(redis, "ECHO HIT:%s:values:%i", key,nVals);
        freeReplyObject(traceRep);
#endif
		for(int i=0; i<nVals; ++i)
		{
			SortableRedisEntry_s entry;
			char *val = (reply->element)[i]->str;
			unpackValue(val, entry.w, entry.f, entry.g);
#ifdef TRACE
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:w:%e:%e:%e:%e:%e:%e:%e", i, entry.w[0], entry.w[1], entry.w[2], entry.w[3], entry.w[4], entry.w[5], entry.w[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:f:%e:%e:%e:%e:%e:%e:%e", i, entry.f[0], entry.f[1], entry.f[2], entry.f[3], entry.f[4], entry.f[5], entry.f[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:g:%e:%e:%e:%e:%e:%e:%e", i, entry.g[0], entry.g[1], entry.g[2], entry.g[3], entry.g[4], entry.g[5], entry.g[6]);
        freeReplyObject(traceRep);
#endif
			entry.distance = 0.0;
			//Compute distance
			for(int j = 0; j < 7; j++)
			{
				entry.distance += fabs(w0[j] - entry.w[j]);
			}
			sortSet.push_back(entry);
		}
		freeReplyObject(reply);
	}

	//Sort the Set/List
	sortSet.sort(SortableRedisEntry_c);

	//Write first bucketSize values from the set to the vectors
	int cnt = 0;
	for(std::list<SortableRedisEntry_s>::iterator iter = sortSet.begin(); iter != sortSet.end(); iter++)
	{
		//Write value to outputs
		double * writeW = new double[7];
		double * writeF = new double[7];
		double * writeG = new double[7];
		for(int i = 0; i < 7; i++)
		{
			writeW[i] = iter->w[i];
			writeF[i] = iter->f[i];
			writeG[i] = iter->g[i];
		}
		wOut->push_back(writeW);
		fOut->push_back(writeF);
		gOut->push_back(writeG);
		//Increment count
		cnt++;
		//Break if we are done
		if(cnt >= bucketSize)
		{
            break;
		}

	}
  delete[] key;
	return;
}


void getCachedSortedSubBucketNearZero(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut, double zeroThresh, std::map<std::string, std::vector<char *> > *dbCache)
{
	std::list<char *> keyList;
	std::vector<double> keyVals[7];
	//Get all values for key
	for(int i = 0; i < 7; i++)
	{
		if( fabs(w0[i]) < zeroThresh && fabs(w0[i] != 0.0))
		{
			(keyVals[i]).push_back(-1.0 * w0[i]);
		}
		keyVals[i].push_back(w0[i]);
	}
	//Funny huge loop time
	for(int i0 = 0; i0 < keyVals[0].size(); i0++)
	{
		for(int i1 = 0; i1 < keyVals[1].size(); i1++)
		{
			for(int i2 = 0; i2 < keyVals[2].size(); i2++)
			{
				for(int i3 = 0; i3 < keyVals[3].size(); i3++)
				{
					for(int i4 = 0; i4 < keyVals[4].size(); i4++)
					{
						for(int i5 = 0; i5 < keyVals[5].size(); i5++)
						{
							for(int i6 = 0; i6 < keyVals[6].size(); i6++)
							{
								double fakeW[7] = { keyVals[0][i0], keyVals[1][i1], keyVals[2][i2], keyVals[3][i3], keyVals[4][i4], keyVals[5][i5], keyVals[6][i6] };
	              char * key = new char[1024];
								buildKey(key, fakeW, tag, keyDigits);
								keyList.push_back(key);
							}
						}
					}
				}
			}
		}
	}

	std::list<SortableRedisEntry_s> sortSet;
	
	//Iterate over keys, and use them
	for(std::list<char *>::iterator iter = keyList.begin(); iter != keyList.end(); iter++)
	{
		//See if we already have this in our map
		std::string mapKey(*iter);
		std::map<std::string, std::vector<char *> >::iterator mapIt = dbCache->find(mapKey);
		if(mapIt != dbCache->end() && !dbCache->empty())
		{
			//We have it
			for(int q = 0; q < mapIt->second.size(); q++)
			{
				//Process point
				SortableRedisEntry_s entry;
				unpackValue(mapIt->second[q], entry.w, entry.f, entry.g);
				entry.distance = 0.0;
				//Compute distance
				for(int j = 0; j < 7; j++)
				{
					entry.distance += fabs(w0[j] - entry.w[j]);
				}
				sortSet.push_back(entry);
			}
		}
		else
		{
			//We did not have it, get it normally
			//Get values
			redisReply *reply;
			reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", *iter);
			assert(reply->type == REDIS_REPLY_ARRAY);
		
			//Write values to set to sort
			int nVals = reply->elements;
			if(nVals == 0)
			{
				freeReplyObject(reply);
#ifdef TRACE
                redisReply * traceRep;
                traceRep = (redisReply *)redisCommand(redis, "ECHO MISS:%s", *iter);
                freeReplyObject(traceRep);
#endif
			}
			else
			{
#ifdef TRACE
                redisReply * traceRep;
                traceRep = (redisReply *)redisCommand(redis, "ECHO HIT:%s:values:%i", *iter,nVals);
                freeReplyObject(traceRep);
#endif
				for(int i=0; i<nVals; ++i)
				{
					SortableRedisEntry_s entry;
					char *val = (reply->element)[i]->str;
					unpackValue(val, entry.w, entry.f, entry.g);
#ifdef TRACE
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:w:%e:%e:%e:%e:%e:%e:%e", i, entry.w[0], entry.w[1], entry.w[2], entry.w[3], entry.w[4], entry.w[5], entry.w[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:f:%e:%e:%e:%e:%e:%e:%e", i, entry.f[0], entry.f[1], entry.f[2], entry.f[3], entry.f[4], entry.f[5], entry.f[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:g:%e:%e:%e:%e:%e:%e:%e", i, entry.g[0], entry.g[1], entry.g[2], entry.g[3], entry.g[4], entry.g[5], entry.g[6]);
        freeReplyObject(traceRep);
#endif
					entry.distance = 0.0;
					//Compute distance
					for(int j = 0; j < 7; j++)
					{
						entry.distance += fabs(w0[j] - entry.w[j]);
					}
					sortSet.push_back(entry);
					//Add entry to map
					char * mapVal = (char *) new double[21];
					memcpy(mapVal, val, sizeof(double) * 21);
					(*dbCache)[mapKey].push_back(mapVal);
				}
				freeReplyObject(reply);
			}
		}

	}

	//Sort the Set/List
	sortSet.sort(SortableRedisEntry_c);

	//Write first bucketSize values from the set to the vectors
	int cnt = 0;
	for(std::vector<double *>::iterator iter = wOut->begin(); iter != wOut->end(); iter++)
	{
    delete[] *iter;
  }
	for(std::vector<double *>::iterator iter = fOut->begin(); iter != fOut->end(); iter++)
	{
    delete[] *iter;
  }
	for(std::vector<double *>::iterator iter = gOut->begin(); iter != gOut->end(); iter++)
	{
    delete[] *iter;
  }
  wOut->clear();
  fOut->clear();
  gOut->clear();
	for(std::list<SortableRedisEntry_s>::iterator iter = sortSet.begin(); iter != sortSet.end(); iter++)
	{
		//Write value to outputs
		double * writeW = new double[7];
		double * writeF = new double[7];
		double * writeG = new double[7];
		for(int i = 0; i < 7; i++)
		{
			writeW[i] = iter->w[i];
			writeF[i] = iter->f[i];
			writeG[i] = iter->g[i];
		}
		wOut->push_back(writeW);
		fOut->push_back(writeF);
		gOut->push_back(writeG);
		//Increment count
		cnt++;
		//Break if we are done
		if(cnt >= bucketSize)
		{
       break;
		}
  }
	for(std::list<char *>::iterator iter = keyList.begin(); iter != keyList.end(); iter++)
	{
    delete[] *iter;
  }
  keyList.clear();
  sortSet.clear();
	return;

}


void getSortedSubBucketNearZero(double w0[7], char * tag, redisContext *redis, int keyDigits, int bucketSize, std::vector<double *> * wOut, std::vector<double *> * fOut, std::vector<double *> * gOut, double zeroThresh)
{
	std::list<char *> keyList;
	std::vector<double> keyVals[7];
	//Get all values for key
	for(int i = 0; i < 7; i++)
	{
		if( fabs(w0[i]) < zeroThresh && fabs(w0[i] != 0.0))
		{
			(keyVals[i]).push_back(-1.0 * w0[i]);
		}
		keyVals[i].push_back(w0[i]);
	}
	//Funny huge loop time
	for(int i0 = 0; i0 < keyVals[0].size(); i0++)
	{
		for(int i1 = 0; i1 < keyVals[1].size(); i1++)
		{
			for(int i2 = 0; i2 < keyVals[2].size(); i2++)
			{
				for(int i3 = 0; i3 < keyVals[3].size(); i3++)
				{
					for(int i4 = 0; i4 < keyVals[4].size(); i4++)
					{
						for(int i5 = 0; i5 < keyVals[5].size(); i5++)
						{
							for(int i6 = 0; i6 < keyVals[6].size(); i6++)
							{
								double fakeW[7] = { keyVals[0][i0], keyVals[1][i1], keyVals[2][i2], keyVals[3][i3], keyVals[4][i4], keyVals[5][i5], keyVals[6][i6] };
	              char * key = new char[1024];
								buildKey(key, fakeW, tag, keyDigits);
								keyList.push_back(key);
							}
						}
					}
				}
			}
		}
	}

    //printf("size %i\n", int(keyList.size()));
	std::list<SortableRedisEntry_s> sortSet;
	//Iterate over keys, and use them
	for(std::list<char *>::iterator iter = keyList.begin(); iter != keyList.end(); iter++)
	{
		//Get values
		redisReply *reply;
		reply = (redisReply *)redisCommand(redis, "SMEMBERS %s", *iter);
		assert(reply->type == REDIS_REPLY_ARRAY);
	
		//Write values to set to sort
		int nVals = reply->elements;
		if(nVals == 0)
		{
			freeReplyObject(reply);
#ifdef TRACE
            redisReply * traceRep;
            traceRep = (redisReply *)redisCommand(redis, "ECHO MISS:%s", *iter);
            freeReplyObject(traceRep);
#endif
		}
		else
		{
#ifdef TRACE
            redisReply * traceRep;
            traceRep = (redisReply *)redisCommand(redis, "ECHO HIT:%s:values:%i", *iter,nVals);
            freeReplyObject(traceRep);
#endif
			for(int i=0; i<nVals; ++i)
			{
				SortableRedisEntry_s entry;
				char *val = (reply->element)[i]->str;
				unpackValue(val, entry.w, entry.f, entry.g);
#ifdef TRACE
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:w:%e:%e:%e:%e:%e:%e:%e", i, entry.w[0], entry.w[1], entry.w[2], entry.w[3], entry.w[4], entry.w[5], entry.w[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:f:%e:%e:%e:%e:%e:%e:%e", i, entry.f[0], entry.f[1], entry.f[2], entry.f[3], entry.f[4], entry.f[5], entry.f[6]);
        freeReplyObject(traceRep);
        traceRep = (redisReply *)redisCommand(redis, "ECHO FETCHED:DATA:%i:g:%e:%e:%e:%e:%e:%e:%e", i, entry.g[0], entry.g[1], entry.g[2], entry.g[3], entry.g[4], entry.g[5], entry.g[6]);
        freeReplyObject(traceRep);
#endif
				entry.distance = 0.0;
				//Compute distance
				for(int j = 0; j < 7; j++)
				{
					entry.distance += fabs(w0[j] - entry.w[j]);
				}
				sortSet.push_back(entry);
			}
			freeReplyObject(reply);
		}
	}

	//Sort the Set/List
	sortSet.sort(SortableRedisEntry_c);

	//Write first bucketSize values from the set to the vectors
  //reset stuff
	int cnt = 0;
  wOut->clear();
  fOut->clear();
  gOut->clear();
  //printf("sortsetsize %i\n", int(sortSet.size()));

	for(std::list<SortableRedisEntry_s>::iterator iter = sortSet.begin(); iter != sortSet.end(); iter++)
	{
		//Write value to outputs
		double * writeW = new double[7];
		double * writeF = new double[7];
		double * writeG = new double[7];
		for(int i = 0; i < 7; i++)
		{
			writeW[i] = iter->w[i];
			writeF[i] = iter->f[i];
			writeG[i] = iter->g[i];
//    printf("iter->w[%i] %.16f\n", i,iter->w[i]);
		}
		wOut->push_back(writeW);
		fOut->push_back(writeF);
		gOut->push_back(writeG);
		//Increment count
		cnt++;
		//Break if we are done
		if(cnt >= bucketSize)
		{
            break;
		}
  }
  //  printf("list finished\n");
	for(std::list<char *>::iterator iter = keyList.begin(); iter != keyList.end(); iter++)
	{
    delete[] *iter;
  }
  keyList.clear();
  sortSet.clear();
	return;
}

std::set<char *> * rebuildDB(std::set<char *> * keyDict, char * tag, redisContext *redis, int newKeyDigits)
{
	std::set<char *> *retSet = new std::set<char *>;
	//Iterate over the old dictionary
	for(std::set<char *>::iterator iter = retSet->begin(); iter != retSet->end(); iter++)
	{
		//Grab bucket corresponding to key
		std::vector<double *> w;
		std::vector<double *> f;
		std::vector<double *> g;
		getBucket(*iter, redis, &w, &f, &g);
		//Re-insert bucket
		for(int i = 0; i < w.size(); i++)
		{
			putData(w[i], f[i], g[i], tag, redis, newKeyDigits);
	    char * newKey = new char[1024];
			buildKey(newKey, w[i], tag, newKeyDigits);
			retSet->insert(newKey);
      delete[] newKey;
		}
		//Delete bucket
		redisReply *reply;
		reply = (redisReply *) redisCommand(redis, "DEL %s", *iter);
		freeReplyObject(reply);
    //FIXME free ok or not?
    freeClear(w);
    freeClear(f);
    freeClear(g);

	}

	return retSet;
}

