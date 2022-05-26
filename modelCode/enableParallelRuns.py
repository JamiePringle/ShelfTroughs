#this module manages a database that allows the parallel running of a
#series of a model over many parts of a parameter space.  The idea is
#that multiple processes suck parameters out of the database, and
#store the results in the database. All of these parameters can run on
#seperate processors.


import sqlite3 as sql
import pickle
import time

#what cPickle protocol to use 
PickleProt=-1

def initializeDB(filename,Plist):
    '''
    If the database in filename exists, this routine does nothing but
    open the database and return the connection to the database.  If
    the database does not exist, it creates the database.

    Plist is a list of tuples, carrying the parameters for the model.
    If the database exists, check to see if all values of Plist exist
    in the database.  If they do not, add them.

    the table is called paramSets, and the consists of the tuple:Pval,
    the value[s] of the parameter to be varied; the boolean isTaken,
    which is False if no one has started to work on this parameter;
    the boolean isDone, which is False until a model has processed the
    data; and the blob:Answer, where Answer is a binary blob that
    contains the dictionary answer that has all the answers :-) in it!
    '''
    #now open connection to database,
    #THERE IS A BUG LURKING HERE IF TWO PROCESSES TRY TO CREATE THE DB AT ONCE!
    conn = sql.connect(filename)
    cur = conn.cursor()

    #check if table exists.  if not, create it, and populate Pval with the values in Plist
    isDone=False
    while not(isDone):
        try:
            cur.execute("pragma table_info(paramSets)")
            isDone=True
        except:
            print('Waiting in pragma table')
            time.sleep(2)
    tableInfo=cur.fetchall()
    print('Table info is' , tableInfo)
    if tableInfo==[]:
        print('Creating table')
        cur.execute( '''CREATE TABLE paramSets (Pval blob, isTaken boolean, isDone boolean, Answers blob)''' )
        for p in Plist:
            cur.execute('INSERT INTO paramSets VALUES (?,?,?,?)', 
                        (sql.Binary(pickle.dumps(p,protocol=PickleProt)),False,False,None))
        conn.commit()
    else:
        print('Table exists already - checking to see if all parameters exist, after 1 second pause')
        #pause is just to make sure someone is already creating the table and filling it in
        time.sleep(2)

        #we want this transaction to be atomic, because it is
        #quite likely that many codes will be trying to do the
        #same thing at the same time; begin exclusive may fail, so
        #be prepared to repeat till succecess
        haveDataBase=False
        while not(haveDataBase):
            try:
                cur.execute('BEGIN EXCLUSIVE')
                haveDataBase=True
            except:
                print('Did not get database... sleeping')
                time.sleep(10.0)

        for p in Plist:
            #only insert if does not exist
            p_pickled=pickle.dumps(p,protocol=PickleProt)                    
            cur.execute('select Pval FROM ParamSets WHERE Pval=?',(sql.Binary(p_pickled),))
            rows=cur.fetchall()
            assert len(rows)<2,'WHAT?!?!? Multiple copies of same parameter in table!?!'
            if len(rows)==0:
                cur.execute('INSERT INTO paramSets VALUES (?,?,?,?)', (sql.Binary(p_pickled),False,False,None))
                print('Inserted')
            else:
                print('Parameter ',' Already exists in table')
                #The two print statements aboe had p in print but I deleted 
        conn.commit()

    return conn

def getNextParam(conn):
    '''
    Takes the connection arguement to the parameter database as an
    arguement, and returns the tuple of the next set of parameters to
    be run.  Marks the run as taken in the database.

    If no more are to be run, returns [].
    '''
    cur=conn.cursor()
    cur.execute('BEGIN EXCLUSIVE')
    cur.execute('SELECT Pval FROM paramSets WHERE isTaken=? Limit 1', (False,))
    rows=cur.fetchall()
    if len(rows)==0:
        pSet=[]
        print('All done, no more parameters left to run')
    else:
        #get first free parameter, and then mark it as taken
        pBlob=rows[0][0]
        pSet=pickle.loads(pBlob)
        cur.execute('update paramSets set isTaken=? where Pval=?',(True,sql.Binary(pBlob))) 

    conn.commit()								
    #print('Working on Pval=',pSet)
    return pSet

def cleanZombies(filename):
    '''warning -- this is a dangerous command.  It is meant to be used
    when something has interrupted runs, and isTaken is true, but the
    run will never be finished.  If finds all places wher isTaken is
    True and isDone is false, and sets isTaken to False.  This should
    allow one to rerun the code to pick up where something went wrong.
    '''
    #now open connection to database,
    conn = sql.connect(filename)
    cur = conn.cursor()

    #now clean things up
    cur.execute('update paramSets set isTaken=? where isTaken=? AND isDone=?',(False,True,False)) 
    conn.commit()
    conn.close()


def finishRun(conn,pSet,answer):
    '''
    Takes the connection to the parameter database, the parameterset and a dictionary
    containing the results of the model, and saves them to the
    database and marks the run as complete.
    '''
    cur=conn.cursor()
    s=pickle.dumps(pSet,protocol=PickleProt)

    #again, if lots of runs are hammering on the data set, this might fail, so try again
    writeWorked=False
    while not(writeWorked):
        try:
            cur.execute('update paramSets set Answers=?,isDone=? where Pval=?',
                        (sql.Binary(pickle.dumps(answer,protocol=PickleProt)),True,sql.Binary(s)))
            conn.commit()
            writeWorked=True
        except:
            print('Write failed, sleep 1 second and try again')
            time.sleep(1.0)

    
def getAllData(filename):
    '''
    Given the filename, return a list of tuples formed form the parameter set
    and the results from the model run, in the form
    (pSet,answers), where pSet is a tuple of parameter values, and
    answers is a dictionary with model results.

    Only returns parameter values and results for parameter sets that
    the models have finished with.

    If the data set is empty, or if the data file does not exist, return empty list. 
    '''
    #open connection to database,
    conn = sql.connect(filename)
    cur = conn.cursor()

    #check if table exists.  if not, complain
    cur.execute("pragma table_info(paramSets)")
    tableInfo=cur.fetchall()
    #print 'Table info is' , tableInfo
    if tableInfo==[]:
        noData=True
        print('There is no file, or no table in the file', filename)
    else:
        #get all the data
        noData=False
        cur.execute('select * from paramSets')
        rows=cur.fetchall()

    #Now loop over all the returned data
    toReturn=[]
    if not(noData):
        for r in rows:
            if r[2]:
                pVal=pickle.loads(r[0])
                Answers=pickle.loads(r[3])
                toReturn.append((pVal,Answers))

        conn.close()

        #now sort data by first element in pVal tuple
        def sortby(p1):
            return p1[0][0]
        toReturn=sorted(toReturn, key = sortby)

    return toReturn

def getSomeData(filename,DataList):
    '''
    Given the filename, return a list of tuples formed form the parameter set
    and the results from the model run, in the form
    (pSet,answers), where pSet is a tuple of parameter values, and
    answers is a dictionary with model results.

    ONLY includes parameters given in DataList, a list of strings

    Only returns parameter values and results for parameter sets that
    the models have finished with.

    If the data set is empty, or if the data file does not exist, return empty list. 
    '''
    #open connection to database,
    conn = sql.connect(filename)
    cur = conn.cursor()

    #check if table exists.  if not, complain
    cur.execute("pragma table_info(paramSets)")
    tableInfo=cur.fetchall()
    print('Table info is' , tableInfo)
    if tableInfo==[]:
        noData=True
        print('There is no file, or no table in the file', filename)
    else:
        #get all the data
        noData=False
        cur.execute('select * from paramSets')
        rows=cur.fetchall()

    #Now loop over all the returned data
    toReturn=[]
    if not(noData):
        for r in rows:
            if r[2]:
                pVal=pickle.loads(r[0])
                Answers_in=pickle.loads(r[3])
                Answers={}
                for ans in DataList:
                    Answers[ans]=Answers_in[ans]
                toReturn.append((pVal,Answers))

        conn.close()

        #now sort data by first element in pVal tuple
        def sortby(p1):
            return p1[0][0]
        toReturn=sorted(toReturn, key = sortby)

    return toReturn

#now lets make a simple test function
if __name__=="__main__":

    fileName='jnk.sqlite' #what file to keep results in. 
    Njobs=20 # number of things to do
    nRuns=4  # number of processes to do it on

    #make parames
    Plist=[(n,) for n in range(Njobs)]

    #initialize paralel runs
    conn=initializeDB(fileName,Plist)
    print('initialized run')
    
    #do data processing
    while True:
        param=getNextParam(conn)
        if param==[]:
            print('ALL DONE WITH ALL parameters')
            break
        else:
            print('Working on param %d'%(param[0],))
            time.sleep(1.0)
            x=param[0]
            answer={'ans':x**2.0}
            finishRun(conn,param,answer)
            print('   Done with job for param %d'%(param[0]))
    
