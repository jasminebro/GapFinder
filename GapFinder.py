import operator,re,linecache

fiveprime=re.compile('pad=(\d+)') #regular expression to pull the start and end bp number for the alu of interest 
threeprime=re.compile('pad=(\d+)') 


PslFile= open('C:\\Users\\jbro262\\Desktop\\SMReferenceAlu\\FUllLenAlu\\SMFullAlu1NoDups.psl',"r") #file with your blat output 
SharedFile=open('C:\\Users\\jbro262\\Desktop\\SMReferenceAlu\\FUllLenAlu\\SMMarRefShared.psl',"a") #file for shared alu
#MarmosetSpecFile=open('C:\\Users\\jbro262\\Desktop\\NWMProject\\PslAnalyzed\\CebusMarSharedVsSM\\SMSpecific_Cebus.psl',"a") #file for Ref. Spec alu
#AotusSpecFile=open('C:\\Users\jbro262\\Desktop\\NWMProject\\PslAnalyzed\\AotusGapsRun2\\AotusSpecificVsMar.psl',"a") #file for Query Spec Alu
AotusInterestFile=open('C:\\Users\\jbro262\\Desktop\\SMReferenceAlu\\FUllLenAlu\\SMSpecificVsMarRef.psl',"a") #Alu of interest,ones I started with 
#AotusBlatAluFile=open('C:\\Users\\jbro262\\Desktop\\NWMProject\\PslAnalyzed\\CebusMarSharedVsSM\\CebusMarBlatAlu_SM.psl',"a")#Alu determined by blat during alignment(approx alu length must make sure it's an alu)
CheckAgainFile=open('C:\\Users\\jbro262\\Desktop\\DuplicatesRemoved\\CategoryFiles\\AotusSMMarVs.Human\\AotusSMMarSharedVsHuman_CheckAgain50.psl',"a") #not sure what these are check again 
Count=open('C:\\Users\\jbro262\\Desktop\\DuplicatesRemoved\\CategoryFiles\\AotusSMMarVs.Human\\CategoryCount.txt',"a") #txt file to keep up with how many of each there are



#PSL COLUMNS --> [0] = Match, [1] = mismatch, [2] = rep. Match, [3] = N's, [4] = Q Gap count, [5] = Qgap bases
# [6] = T Gap count, [7] = T gap bases, [8] = Strand, [9] = Q name, [10] = Q size, [11] = Q Start
#[12] = Q end, [13] = T name, [14] = Tsize, [15] = T start, [16] = T end, [17] = BlockCount, [18] = Block Sizes
# [19] = qStarts, [20] = tStarts  

QueryGapStarts=[]
QueryGapSizes=[]
TargetGapStarts=[]
TargetGapSizes=[]
AluInterestStart=[]
ChrLocations=[]
AluStartLocation=[]
AluEndLocation=[]
SharedList=[]
MarmosetSpecList=[]
AotusSpecList=[] 
AluOfInterest=[]
AluFromBlat=[]
CheckAgainList=[]


for line in PslFile:
    column=line.split('\t')
    left=(fiveprime.findall(column[9])) #the number in list format
    right=(threeprime.findall(column[9]))
    AluInterestLeft = [int(x) for x in left]
    AluInterestStart.append(AluInterestLeft)
    AluInterestRight=[int(x) for x in right]
    QueryBlockStarts=column[19].split(',')
    BlockSizes=column[18].split(',')
    TargetBlockStarts=column[20].split(',')
    ChrLocation=column[13] 
    StartInChr=column[15]
    EndInChr=column[16]
    del QueryBlockStarts[-1]
    del BlockSizes[-1]
    del TargetBlockStarts[-1]
    QueryBlockStarts=list(map(int,QueryBlockStarts))
    TargetBlockStarts=list(map(int,TargetBlockStarts))
    BlockSizes=list(map(int,BlockSizes))
    GapStartQuery=list(map(operator.add,QueryBlockStarts,BlockSizes))
    #print(QueryBlockStarts)
    #print(BlockSizes)
    #print(GapStartQuery)
    GapSizeQuery=[x-y for x,y in zip(QueryBlockStarts[1:],GapStartQuery)]  
    QueryGapStarts.append(GapStartQuery)
    QueryGapSizes.append(GapSizeQuery)                
    GapStartTarget=list(map(operator.add,TargetBlockStarts,BlockSizes))
    GapSizeTarget=[x-y for x,y in zip(TargetBlockStarts[1:],GapStartTarget)] 
    TargetGapStarts.append(GapStartTarget)
    TargetGapSizes.append(GapSizeTarget)
    ChrLocations.append(ChrLocation)
    AluStartLocation.append(StartInChr) 
    AluEndLocation.append(EndInChr) 


"""
These list will hold the line numbers that I need to pull from the file
But I need to keep track of which alu I am looking at by keeping track of the index value in the GapStarts list
Maybe create two lists for AluOfInterest line numbers and NewlyFoundAlu line numbers
i= the line number in the psl file(also the numbering of the Alu start and end numbers you should have enough start numbers for each
alu therefore the number of lines and that list length should be the same 
"""

for i, (x,y,z) in enumerate(zip(TargetGapSizes,QueryGapSizes,QueryGapStarts)): #keeps all 3 lists together as the loop loops
    
    if all (m<50 for m in x) and all(n <50 for n in y): #Shared files sorted here
        SharedList.append(i+1)
        line=linecache.getline('C:\\Users\\jbro262\\Desktop\\SMReferenceAlu\\FUllLenAlu\\SMFullAlu1NoDups.psl',SharedList[-1])
        SharedFile.write("TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+ str(y)+'\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"SharedAluStartInMarmoset:"+ str(AluStartLocation[i])+'\t'+"SharedAluEndInMarmoset:"+str(AluEndLocation[i])+'\n'+line)
    #elif any(330>=o>=260 for o in x) and all(p<260 for p in y): #Ref Spec files sorted here 
        #MarmosetSpecList.append(i+1)
        #line=linecache.getline('C:\\Users\\jbro262\\Desktop\\CebusMarSharedBlatSMFiltered.psl',MarmosetSpecList[-1])
        #MarmosetSpecFile.write("TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+ str(y)+ '\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"AluStartInCebus:"+ str(AluStartLocation[i])+'\t'+"AluEndInCebus:"+str(AluEndLocation[i])+'\n'+line)

    elif all(q<50 for q in x) and any(330>=r>=260 for r in y): #Query Spec files sorted here 
        BlatAluStart=[k for (k,j) in zip(z,y) if 330>=j>=260]
        #AotusSpecList.append(i+1) 
        #line=linecache.getline('C:\\Users\\jbro262\\Desktop\\CebusMarSharedBlatSMFiltered.psl',AotusSpecList[-1])
        #AotusSpecFile.write("TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+ str(y)+ '\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"MissingAluStartInSM:"+ str(AluStartLocation[i])+'\t'+"MissingAluEndInSM:"+str(AluEndLocation[i])+'\n'+line)

        for num in BlatAluStart: #checking location to see if it's the Alu I started with or something else detected during alignments
            if  570 <=num<=630: #30 basepairs on both sides of the Alu of interst 
                AluOfInterest.append(i+1)
                line=linecache.getline('C:\\Users\\jbro262\\Desktop\\SMReferenceAlu\\FUllLenAlu\\SMFullAlu1NoDups.psl',AluOfInterest[-1])
                AotusInterestFile.write("AluOfInterestStartBp:"+str(AluInterestStart[i])+'\t'+"BlatAluStart:"+str(BlatAluStart)+'\t'+"TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+str(y)+ '\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"AotusAluInterestStartInMar:"+str(AluStartLocation[i])+'\t'+"AotusAluInterestEndInSM:"+str(AluEndLocation[i])+'\n'+line)                   
            #else: #alu detected during alignments go here 
                #AluFromBlat.append(i+1)
                #line=linecache.getline('C:\\Users\\jbro262\\Desktop\\CebusMarSharedBlatSMFiltered.psl',AluFromBlat[-1])
                #AotusBlatAluFile.write("AluOfInterestStartBp:"+str(AluInterestStart[i+1])+'\t'+"BlatAluStart:"+str(BlatAluStart)+'\t'+"TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+str(y)+ '\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"AoAluBlatStartInSM:"+ str(AluStartLocation[i])+'\t'+"AotusAluBlatEndInSM:"+str(AluEndLocation[i])+'\n'+line)            

    else: #Check these are something funky is going on 
        CheckAgainList.append(i+1) 
        line=linecache.getline('C:\\Users\\jbro262\\Desktop\\DuplicatesRemoved\\CategoryFiles\\AotusSMMarVs.Human\\AotusMarSMSharedVsHumanClean.psl',CheckAgainList[-1])
        CheckAgainFile.write("TargetGapSizes:"+ str(x)+'\t'+"QueryGapSizes:"+str(y)+ '\t'+"QueryGapStarts:"+str(z)+'\t'+"ChrNum:"+str(ChrLocations[i])+'\t'+"CheckAgainStartInSM:"+ str(AluStartLocation[i])+'\t'+"CheckAgainEndInSM:"+str(AluEndLocation[i])+'\n'+line)

SharedFile.close()
#MarmosetSpecFile.close()
#AotusSpecFile.close()
AotusInterestFile.close()
#AotusBlatAluFile.close() 
CheckAgainFile.close()     
    
Count.write("Number of AotusSMMarHuman Shared Queries:"+str(len(SharedList))+'\n')
#Count.write("Number of Squirrel Monkey Specific Queries:"+ str(len(MarmosetSpecList))+'\n')
#Count.write("Number of AotusMar_SM Specific Queries:"+ str(len(AotusSpecList))+'\n')
Count.write("Number of AotusMarSM Shared Alu :" + str(len(AluOfInterest))+'\n')
#Count.write("Number of Alu from Blat Identified:" + str(len(AluFromBlat))+'\n') 
Count.write("Number of AotusMarSM Shared Alu to Check Again:"+ str(len(CheckAgainList))+'\n') 
Count.close() 
   
print("DONE!!!") 

