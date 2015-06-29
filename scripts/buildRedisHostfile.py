import os
import sys

def genNodeList():
    #Get environmental variable
    nodeString = os.environ['SLURM_JOB_NODELIST']
    #For debug/devel
    #nodeString = "cn[14,3-9,46,52,61,111-113,116-118]"
    #Remove leading 'cn'
    nodeString = nodeString.lstrip('cn')
    nodeList = []
    if (nodeString.startswith('[') == False):
        #Just the one
        nodeList.append(int(nodeString))
    else:
        #Strip the brackets, we no longer need them
        nodeString = nodeString.strip('[]')
        #Tokenize the string
        tokenList = nodeString.rsplit(',')
        #Iterate over tokens
        for token in tokenList:
            #Two cases: Single or range
            if '-' in token:
                inclusiveRange = token.rsplit('-')
                start = int(inclusiveRange[0])
                end = int(inclusiveRange[1]) + 1
                for node in xrange(start, end):
                    nodeList.append(node)
            else:
                nodeList.append(int(token))
    return nodeList

#Attempts to break nodeList into numServers chunks
#Does not guarantee numServers chunks, but will attempt to get close
def assignServers(nodeList, numServers):
    serverPairs = []
    serverList = []
    #See if this is one of the easy case
    #TODO: Do we need the cases anymore now that we have smart-ish chunking?
    if (numServers == 1):
        #Just one server, so return the first
        for node in nodeList:
            serverPairs.append( (node, nodeList[0]) )
        serverList.append(nodeList[0])
    elif ( numServers >= len(nodeList) ):
        #Everyone is a server
        for node in nodeList:
            serverPairs.append( (node, node) )
        serverList = nodeList
    else:
        #We need to assign servers
        #TODO: Make a more even distribution
        #TODO: Guarantee we have numServers chunks
        numNodes = len(nodeList)
        nodesPerServer = numNodes / numServers
        remNodes = numNodes % numServers
        if(remNodes > 0):
            nodesPerServer = nodesPerServer + 1
        chunkList = [nodeList[i:i + nodesPerServer] for i in xrange(0, len(nodeList), nodesPerServer)]
        #chunks[i][0] is a server
        for chunk in chunkList:
            serverList.append(chunk[0])
            for node in chunk:
                serverPairs.append( (node, chunk[0]) )
    return (serverList, serverPairs)

def writeLUTFile(serverPairs):
    #Open file to write
    lFile = open('lutFile', 'w')
    for pair in serverPairs:
        line = "cn" + str(pair[0]) + "\t" + "cn" + str(pair[1]) + "\n"
        lFile.write(line)
    lFile.close()

def writeHostFile(nodeList):
    nFile = open ('nodeFile', 'w')
    for node in nodeList:
        line = "cn" + str(node) + "\n"
        nFile.write(line)
    nFile.close()

def writeServerFile(serverList):
    sFile = open ('serverFile', 'w')
    for node in serverList:
        line = "cn" + str(node) + "\n"
        sFile.write(line)
    sFile.close()

def writeRedisConfigs(serverList):
    pwd = os.getcwd()
    master = "cn" + str(serverList[0])
    for node in serverList:
        nHost = "cn" + str(node)
        #See if we need to prepare the child
        if(nHost != master):
            #Make child dir if needed
            if not os.path.exists(nHost):
                os.mkdir(nHost)
            #Go to child dir
            os.chdir(nHost)
            #Write file
            rFile = open('redis.conf', 'w')
            rFile.write("port 6379\n")
            rFile.write("slaveof " + master + " 6379\n")
            rFile.write("dbfilename dump" + nHost + ".rdb\n")
            rFile.write("slave-read-only no\n")
            rFile.close()
            #Return to parent dir
            os.chdir(pwd)

def which(program):
    for path in os.environ['PATH'].split(os.pathsep):
        path = path.strip('""')
        binary = os.path.join(path, program)
        if os.path.isfile(binary):
            return binary

def writeBashScripts(serverList):
    pwd = os.getcwd()
    master = "cn" + str(serverList[0])
    sFile = open('startRedisServers.sh', 'w')
    eFile = open('endRedisServers.sh', 'w')
    #Boiler plate
    bashLine = "#!/bin/bash\n"
    sFile.write(bashLine)
    eFile.write(bashLine)
    #Do the start script commands
    #Find redis-server
    rServer = which('redis-server')
    for node in serverList:
        nHost = "cn" + str(node)
        if(nHost != master):
            #Go to child dir
            os.chdir(nHost)
            #Get child dir
            kwd = os.getcwd()
            #Write command line
            cLine = "nohup ssh -n " + nHost
            cLine += " \'cd " + kwd + "; " + rServer + " " + kwd + "/redis.conf &\' &\n"
            sFile.write(cLine)
            #Return to parent dir
            os.chdir(pwd)
        else:
            #Master is easy. Also it is probably localhost
            cLine = "nohup ssh -n " + nHost
            cLine += " \'" + rServer + " &' &\n"
            sFile.write(cLine)
    sFile.close()
    #Now write the end script
    killswitch = "pkill -f \"redis-server\""
    for node in serverList:
        nHost = "cn" + str(node)
        cLine = "ssh " + nHost + " \'" + killswitch + "\'\n"
        eFile.write(cLine)
        eFile.write(cLine)
    eFile.close()

def main():
    numServers = 1
    if(len(sys.argv) == 2):
        numServers = int(sys.argv[1])
    nodeList = genNodeList()
    (serverList, serverPairs) = assignServers(nodeList, numServers)
    #Write LUT file for actual program to use
    writeLUTFile(serverPairs)
    #Write nodelist and serverlist if the following scripts aren't used
    writeHostFile(nodeList)
    writeServerFile(serverList)
    #Write redis config files and startup/shutdown scripts
    writeRedisConfigs(serverList)
    writeBashScripts(serverList)


if __name__ == "__main__":
    main()
