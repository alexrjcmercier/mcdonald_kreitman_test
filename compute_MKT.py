#coding:utf-8
#!/bin/env python2.7
#This python script written for python 2.7.x (should still be compatible with older distributions though) generates a McDonald-Kreitman Test table to be used with SnIPRE, with the values of Ps, Pn, Ds, Dn plus a simple MK Test with no correction and a Pn/Ps column for reference to an older test. It handles recusrively as many files as there is in the given directory.
#It requires a genetic code tabulated text file and a multi FASTA file for each gene, with sequences of this gene for every individual as well as its appartenance to a group (@1 or @2), the groups will be indentified as principal or outgroup based on the number of individuals (the ougroup need to be restricted to only 1 individual).
#Arguments must be given in this order "input_dir; genetic_code_file; outputdir"

####################################################################################################################################
#                                              CALLING MODULES AND FETCHING ARGUMENTS                                              #
####################################################################################################################################

import os, sys, itertools
script, inputdir, inputtxtgencode, outputdir = sys.argv

####################################################################################################################################
#                                                       FUNCTIONS DEFINITION                                                       #
####################################################################################################################################

#GetGeneticCode(x) - Fetch the genetic code from a file (given as argument) and store it in a dictionary with codon as key and aminoacid as value.
#Input: the genetic code as a tabulated .txt file with, in this order, the full AA name, the 3 letters AA name, the 1 letter AA name, 1st base of codon, 2nd base, 3rd base (e.g. Phenylalanine\tPhe\tF\tT\tT\tT)
def GetGeneticCode(inputfile) :
    GeneticCode = {}
    STOPCodons = []
    with open(inputfile, "r") as InputGeneticCode :
        Lineread = InputGeneticCode.readline().rstrip().split("\t")
        while Lineread != [''] :
            GeneticCode["".join([str(x) for x in [Lineread[-3], Lineread[-2], Lineread[-1]]])] = Lineread[1]
            Lineread = InputGeneticCode.readline().rstrip().split("\t")
    for key, value in GeneticCode.items() :
        if value == "STP" :
            STOPCodons.append(key)
    return GeneticCode, STOPCodons
    #Output: a dictionary with codons as keys and aminoacids as values and a list of STOP codons.

#ListFiles(x) - List all the files to be analysed, determine the number or studied individuals and which are part of the ingroup or outgroup.
#Input: a directory where files to analyse are located.
def ListFilesAndGroupMembers(inputdir) :
    FilesList = [] #List the files and sort the list alphabetically.
    for File in os.listdir(inputdir) :
        FilesList.append(File)
    FilesList.sort()
    Individuals, NbIndividuals, Ingroup, Outgroup, IngroupKey, OutgroupKey = {}, 0, [], [], '', ''  #List the group members (@1 and @2) from the first open file.
    with open ("".join([inputdir, FilesList[0]]), "r") as InputFile :
        for Line in InputFile :
            if Line[0] == ">" :
                if Line[-2::] not in Individuals :
                    Individuals[Line[-2::]] = []
                Individuals[Line[-2::]].append(Line[1:-4])
    for key, value in Individuals.items() :
        if len(value) == 1 :
            Outgroup, OutgroupKey = value, key
        else :
            Ingroup, IngroupKey = value, key
    Nbindividuals = len(Ingroup) + len(Outgroup)
    return FilesList, Ingroup, Outgroup, NbIndividuals
    #Output: a tuple of variables.

#GenDictOfSequences(x) - Generate a dictionary storing the fasta sequences for all individuals from the open file plus the indication if there is a missing sequence for ingroup and for outgroup members. Warning: resets the dictionary at each call of the function.
#Input: a FASTA-formatted sequences file
def GenDictOfSequences(inputfile) :
    Sequences_dict = {}
    Lenseq = 0
    with open (inputfile, "r") as Inputfile :
        lineread1 = Inputfile.readline().rstrip()
        lineread2 = Inputfile.readline().rstrip()
        if len(lineread2) != 0 :
            Lenseq = len(lineread2)
        while lineread1 != "" :
            if len(lineread2) == Lenseq :
                Sequences_dict[lineread1[1:-3]] = lineread2
            lineread1 = Inputfile.readline().rstrip()
            lineread2 = Inputfile.readline().rstrip()
        if outgroup[0] not in Sequences_dict.keys() :
            Sequences_dict["missing_outgroup"] = 1
        else :
            Sequences_dict["missing_outgroup"] = 0
        Missing_ingroup = 0
        for individual in ingroup :
            if individual not in Sequences_dict.keys() :
                Missing_ingroup += 1
        if Missing_ingroup >= 1 :
            Sequences_dict["missing_ingroup"] = Missing_ingroup
        else :
            Sequences_dict["missing_ingroup"] = 0
    return Sequences_dict, Lenseq
    #Output: a disctionnary storing the sequences with the individual names as keys, and the sequences as values.

#GetInAndOutCodons(x, y, z) - For the given position (codon index), get the codon sequence for the individuals in the ingroup and the outrgoup.
#Input: in this order: the codon index, the sequences dictionary, the outgroup individual (determine the ingroup by elimination).
def GetInAndOutCodons(codon, sequences_dict, outgroup) :
    Codon_in_Dict, Codon_out_Dict = {}, {}
    thirdbase_index = int((codon*3)-1)
    secondbase_index = thirdbase_index-1
    firstbase_index = thirdbase_index-2
    for key, value in sequences_dict.items() :
        if key not in ["missing_ingroup", "missing_outgroup"] and key != outgroup[0] :
            Codon_in_Dict[key] = str(sequences_dict[key][firstbase_index]+sequences_dict[key][secondbase_index]+sequences_dict[key][thirdbase_index])
        if key not in ["missing_ingroup", "missing_outgroup"] and key == outgroup[0] :
            Codon_out_Dict[key] = str(sequences_dict[key][firstbase_index]+sequences_dict[key][secondbase_index]+sequences_dict[key][thirdbase_index])
    return Codon_in_Dict, Codon_out_Dict
    #Output; two dictionaries with the ingroup codons and outgroup codon.    

#PrintWarningMessages(x, y, z) - Function to print a warning or error message depending on the given warning type.
#Input: the file the script is iterating on; the codon the file is iterating on; the type of warning or error raised.
def PrintWarningMessages(file, codon, warningtype) :
    if warningtype == "MD" :
        print("".join([str(x) for x in ["Warning - Missing data for at least one individual in gene file: ", file, ". Removing file from analysis."]]))
    elif warningtype == "UE" :
        print("".join([str(x) for x in ["Error at codon: ", codon, " of gene file: ", file, " - Unmanaged exception."]]))
    elif warningtype == "NP" :
        print("".join([str(x) for x in ["Warning - No valid parcimonious path possible at codon: ", codon, " of gene file: ", file, " (STOP codon not bordered by STOPs on both sides). Removing codon from analysis."]]))
    elif warningtype == "STOP1" :
        print("".join([str(x) for x in ["Warning - STOP codon detected at the last position (codon: ", codon, ") of gene file: ", file, ". STOP codon removed from analysis."]]))
    else :
        print('Error - Unmanaged exception in function "PrintWarningMessages". Ironic, isn\'t it ?')
    #Output: a printed message. this function does not return anything.

#STOPCount(x, y, z) - Count the number of STOP codons in a list of codons.
#Input: list of possible STOP codons; list of codons to count in; 
def STOPCount(STOPcodons, totcodonslist) :
    STOPDetected = 0
    for STOPCodon in STOPcodons :
        if STOPCodon in totcodonslist :
            STOPDetected = 1
    return STOPDetected
    #Output: a flag with a value of 1 if a STOP was detected, 0 if none detected.    

#CountBaseDiff(x) - In a permutation pair, determine the number of bases difference.
#Input: a permutation pair as a list of two codons.    
def CountBaseDiff(permut_pair) :
    BaseDiff = 0
    for pos in range(0, 3) :
        if permut_pair[0][pos] != permut_pair[1][pos] :
            BaseDiff += 1
    return BaseDiff
    #Output: the number of differences as a integer.

#IncrementPsPnDsDn(a, b, c, d, e, f, g, h) - for a permutation pair, increment the global variables given (filePs, Pn, Ds, Dn) by detecting the synonymous or non-synonymous base differences. Needs to be told if these differences are a polymorphisms or divergences (using the flags inP and inD, one must be equal to 1 and the other 0). Also needs to be given the genetic code dictionary.
#Input: in this order, genetic code dictionary; permutation pair; flag for divergence; flag for polymorphism; global variable for number of Ps; Pn; Ds; Dn.
def IncrementPsPnDsDn(gencode_dict, permut_pair, inP, inD, filePs, filePn, fileDs, fileDn) :
    ErrorMessage = 'Error - Umanaged exception in function "IncrementPsPnDsDn".'
    if inP == 1 and inD == 0 :
        if gencode_dict[permut_pair[0]] == gencode_dict[permut_pair[1]] :
            filePs += 1
        elif gencode_dict[permut_pair[0]] != gencode_dict[permut_pair[1]] :
            filePn += 1
        else :
            print(ErrorMessage)
    elif inP == 0 and inD == 1 :
        if gencode_dict[permut_pair[0]] == gencode_dict[permut_pair[1]] :
            fileDs += 1
        elif gencode_dict[permut_pair[0]] != gencode_dict[permut_pair[1]] :
            fileDn += 1
        else :
            print(ErrorMessage)
    else :
        print(ErrorMessage)
    return filePs, filePn, fileDs, fileDn
    #Output: the incremented global variables; may need to be reassigned (e.g. globvar += (function(initialglobvar))).

#GetBestIntermediaryCombination(x) - Selects the best possible intermediary codon(s) when there is more than one base difference between the members of a given pair of codons. Also handles the filtering of STOP codons.
#Warning - Uses the function PrintWarningMessages in its core. Said function must be declared beforehand.
#Input: a pair of codons
def GetBestIntermediaryCodons(codon_pair, STOPcodons) :
    CombinationList = [[], [], []] #Initiate a list of lists which will contain the different possible nucleotides at the three positions of the codon.
    for pos in range(0, 3) : #Python indexes! pos will take values of 0, 1 and 2
        for PairMember in codon_pair :
            if PairMember[pos] not in CombinationList[pos] :
                CombinationList[pos].append(PairMember[pos])
        IntermediaryCodons = list(itertools.product(*CombinationList)) #itertools.product generates all the possible combination of a set of lists. In this case, all the possible codons with the lists of bases at position1, 2 and 3.
        for PairMember in codon_pair :
            if ((PairMember[0], PairMember[1], PairMember[2])) in IntermediaryCodons : #Remove the intermediary codons identical to the two initial codons of the studied pair.
                ToRM = ((PairMember[0], PairMember[1], PairMember[2]))
                IntermediaryCodons.remove(ToRM)
    ToRM = []
    for IntermediaryCodon in IntermediaryCodons : #Eliminating intermediary STOP codons which are not bordered by initial STOP codons.
        LeftBorder, RightBorder = 0, 0
        IntCodon = "".join(IntermediaryCodon)
        if IntCodon in STOPcodons :
            if codon_pair[0] in STOPcodons : #Checking if one of the initial codons is a stop codon
                LeftBorder = 1
            if codon_pair[-1] in STOPcodons :
                RightBorder = 1
            if LeftBorder == 1 and RightBorder == 1 : #If both initial codons are stop codons, the intermediary stop codon is valid.
                print("".join([str(x) for x in ["Warning - Intermediary STOP codon found correctly bordered by initial STOP codons at codon: ", codon, " in gene file: '", file, "'. Kept in list of possible intermediary codons."]]))
            elif (LeftBorder == 1 and RightBorder == 0) or (LeftBorder == 0 and RightBorder == 1) : #If only halfly bordered (or not bordered at all) by a stop codon, invalid intermediary stop codon.
                ToRM.append(IntermediaryCodon)
                print("".join([str(x) for x in ["Warning - Intermediary STOP codon found not bordered on both sides by initial STOP codons at codon: ", codon, " in gene file: '", file, "'. Deleted from list of possible intermediary codons."]]))
            elif (LeftBorder == 0 and RightBorder == 0) :
                ToRM.append(IntermediaryCodon)
                print("".join([str(x) for x in ["Warning - Intermediary STOP codon found not bordered on any side by initial STOP codons at codon: ", codon, " in gene file: '", file, "'. Deleted from list of possible intermediary codons."]]))
            else :
                print('Error - Unmanaged exception in function "GetBestIntermediaryCodons".')
    for CodonToRM in ToRM : #Removing the invalid STOP codons.
        IntermediaryCodons.remove(CodonToRM)
    if len(IntermediaryCodons) >= 1 : #If valid intermediary codons (there is cases with no valid intermediary codons, in which case the codon is removed from analysis)
        IntermediaryPermutations = []
        IntermediaryBorderedPermutations = []
        for i in itertools.permutations(IntermediaryCodons, basediff-1) : #Creating all the possibles orders of intermediary codons (i.e. mutation apparition order)
            IntermediaryPermutations.append(i)
        for i in IntermediaryPermutations : #Adding the borders (i.e. initial two codons of the pair) to all permutations.
            BufferList = []
            BufferList.append(codon_pair[0])
            for j in i :
                BufferList.append("".join(j))
            BufferList.append(codon_pair[1])
            IntermediaryBorderedPermutations.append(tuple(BufferList))
        IntermediaryScoredCombinations = []
        for IntermediaryPermutation in IntermediaryBorderedPermutations : #Scoring the permutations using (i) the number of nucleotide differences and (ii) the number of aminoacind changes.
            IntermediaryTotalBasesDifference = 0
            for i in range(1, len(IntermediaryPermutation)) : #(i)
                IntermediaryNumberBaseDifference = 0
                for j in range(0,3) :
                    if IntermediaryPermutation[i-1][j] != IntermediaryPermutation[i][j] :
                        IntermediaryNumberBaseDifference += 1
                IntermediaryTotalBasesDifference += IntermediaryNumberBaseDifference
            IntermediaryAAList = []
            for IntermediaryAA in IntermediaryPermutation : #(ii)
                IntermediaryAAList.append(gencode_dict[IntermediaryAA])
                NumberOfAASwaps = 0
                for Swap in range(1, len(IntermediaryAAList)) :
                    if IntermediaryAAList[Swap] != IntermediaryAAList[Swap-1] :
                        NumberOfAASwaps += 1
            IntermediaryScoredCombinations.append([IntermediaryTotalBasesDifference, NumberOfAASwaps, IntermediaryPermutation])
        BestIntermediaryCombination = sorted(IntermediaryScoredCombinations)[0][-1] #Sorting by number of number of bases differences then by number of aminoacids changes, then getting the first found best intermediary combination.
    else :
        BestIntermediaryCombination = None
    return BestIntermediaryCombination
    #Output: a tuple of codons with the initial pair as external borders.

#GetPairsInPermut(x) - Isolate all the pairs of codons in the given permutation.
#Input: a permutation of codons.
def GetPairsInPermut(permutation) :
    if permutation != None :
        PermutationPairs = []
        for i in range(1, len(permutation)) :
            PermutationPairs.append([permutation[i-1], permutation[i]])
    else :
        PermutationPairs = None
    return PermutationPairs
    #Output: a list containing pairs as lists of two elements.    

#GetScoredPermutations(x, y) - From a list of codons, list all possible permutations (intermediate steps between these codons) and assign score to it based on the number of nucleotide changes.
#Input: a list of codons.
def GetScoredPermutations(totcodons) :
    CodonPermutations, PermutationsScores, NbTotCodons = [], {}, len(totcodons)
    for Permutation in itertools.permutations(totcodons, NbTotCodons) : #Listing all possible permutations
        CodonPermutations.append(Permutation)
    for Permutation in CodonPermutations : #Counting the number of changes per permutation and building a disctionary with scores assigned to permutations
        TotDiff = 0
        for i in range(1, len(Permutation)) :
            NbDiff = 0
            for j in range(0,3) :
                if Permutation[i-1][j] != Permutation[i][j] :
                    NbDiff += 1
            TotDiff += NbDiff
        PermutationsScores[Permutation] = TotDiff
    return PermutationsScores
    #Output: a dictionary with permutations as keys, and scores as values. The lower the score, the more parcimonious the permutation.

#GetMostParcimoniousPermutation(x, y) - Chose the most parcominious permutation from a dictionary of scored permutations.
#Input: the dictionary of scored mutations, and a list with the outgroup codons (used for parcimony).
def GetMostParcimoniousPermutation(permutations_scores, codon_out_dict) :
    ParcimonyKeys = []
    for key, value in permutations_scores.items() : #Extracting the permutations with the lowest score (lowest number of nucleotide changes).
        if value == min(permutations_scores.values()) :
            ParcimonyKeys.append(key)
    ChosenPermut = ParcimonyKeys[0] #Chosing a preliminary, most parcimonious permutation (first found).
    OutCodonAtEnd = 0
    for key in ParcimonyKeys : #Check for a same-score permutation with outgroup codon at the end. Keeps first found permutation if impossible. Could implement a selection influenced by preferential codons in the gemone.
        if codon_out_dict.values()[0] in key[-1] and OutCodonAtEnd == 0 :
            ChosenPermut = key
            OutCodonAtEnd = 1
    return ChosenPermut, OutCodonAtEnd
    #Output: the chosen most parcimonious permutation as an ordered list of codons, and a flag indication if the ougroup codon is at one end of the permutation (1), or if it isn't (0).

#CalculateNeutralityIndex(w, x, y, z) - Calculate the Neutrality Index as in a basic McDonald-Kreitman Test.
#Input: four variables containint the values of ps, pn, ds and dn in this order.
def CalculateNeutralityIndex(ps, pn, ds, dn) :
    NIBaseElements = []
    for x in [ps, pn, ds, dn] :
        if x < 1 :
            x = 1
        NIBaseElements.append(float(x))
    NI = ((NIBaseElements[1]/NIBaseElements[0])/(NIBaseElements[3]/NIBaseElements[2]))
    return NI
    #Output: the value of Neutrality Index.

#CalculatePnPs(x, y) - Calculate the Pn/Ps ratio.
#Input: two variables containing the values of ps and pn, in this order.
def CalculatePnPs(ps, pn) :
    PnPsElements = []
    for x in [ps, pn] :
        if x < 1 :
            x = 1
        PnPsElements.append(float(x))
    PnPs = ((PnPsElements[1]/PnPsElements[0]))
    return PnPs
    #Output: the value of Pn/Ps ratio.

####################################################################################################################################
#                                                           CORE PROGRAM                                                           #
####################################################################################################################################

gencode_dict, STOPcodons = GetGeneticCode(inputtxtgencode) #Fetching the genetic code from file.
basecode = ["A", "T", "C", "G"]
filelist, ingroup, outgroup, nbindividuals = ListFilesAndGroupMembers(inputdir) #Building a list of files to analyse and determining the number of individuals studied.
with open("".join([outputdir, inputdir.split("/")[-2], "_MKTable.txt"]), "w") as outputfile :
    nbingroups, nboutgroups = len(ingroup), len(outgroup)
    outputlist = ["\t".join(['"geneID"', '"PS"', '"PR"', '"FS"', '"FR"', '"nbstop"', '"NI"','"Pn/Ps"' , '"Tsil"', '"Trepl"', '"Lentgh"', '"nout"', '"npop"', "\n"])] #Initializing the output list with the SnIPRE-formatted header.
    for file in filelist :
        sequences_dict, lenseq = GenDictOfSequences("".join([inputdir, file]))
        filePs, filePn, fileDs, fileDn, fileSTOP, filetotSs, filetotNSs = 0, 0, 0, 0, 0, 0, 0
        totSs, totNSs = 0, 0 #Initializing the total numbers of synonymous and non-synonymous sites.
        if sequences_dict["missing_ingroup"] <= 0 and sequences_dict["missing_outgroup"] <= 0 :
            for codon in range(1, ((lenseq/3)+1)) : #Codon starts at 1.
                gapincodon = 0
                codon_in_dict, codon_out_dict = GetInAndOutCodons(codon, sequences_dict, outgroup)
                incodonssetlist, outcodonssetlist = list(set(codon_in_dict.values())), list(set(codon_out_dict.values()))
                totcodonslist = incodonssetlist+outcodonssetlist
                if "-" in codon_in_dict.values() or "-" in codon_out_dict.values() : #Rejecting codons with one or more gaps.
                    gapincodon = 1
                if gapincodon == 0 :
                    codonSs, codonNSs = 0, 0
                    nbSspos, nbNSspos = 0, 0
                    invarcodon = "".join(list(set(codon_in_dict.values()))[0])
                    invarAA = gencode_dict[invarcodon]
                    for basepos in range(0, 3) :
                        templist = [list(x) for x in [invarcodon, invarcodon, invarcodon]]
                        templistindex = 0
                        tempAAlist = []
                        for base in basecode :
                            if base != invarcodon[basepos] :
                                templist[templistindex][basepos] = base
                                templist[templistindex] = "".join(templist[templistindex])
                                templistindex += 1
                        for tempcodon in templist :
                            tempAAlist.append(gencode_dict[tempcodon])
                        for AA in tempAAlist :
                            if AA == invarAA :
                                nbSspos += 1
                            else :
                                nbNSspos += 1
                        codonSs = (float(nbSspos)/3.0)
                        codonNSs = (float(nbNSspos)/3.0)
                    totSs += codonSs
                    totNSs += codonNSs
                if gapincodon == 1 :
                    print("".join([str(x) for x in ["Warning - Gap detected in codon ", codon, " of gene file ", file, ". Codon removed from analysis."]]))
                    totSs += 1.0
                    totNSs += 2.0
                elif len(list(set(codon_in_dict.values()+codon_out_dict.values()))) == 1 : #Using this test only to eliminate the unvariant sites for the rest of the analysis.
                    pass
                elif codon == ((lenseq/3)+1) : #Eliminating STOP codons at end of sequence.
                    for lastcodon in totcodonlist :
                        if lastcodon in STOPcodons :
                            PrintWarningMessages(file, codon, "STOP1")
                elif len(list(set(codon_in_dict.values()+codon_out_dict.values()))) == 2 : #Checking for polymorphisms and divergences: simpler case where only two codons are possible in all sequences (ingroup+outgroup).
                    STOPdetected = STOPCount(STOPcodons, totcodonslist)
                    fileSTOP += STOPdetected
                    if len(set(codon_in_dict.values())) == 2 : #Polymorphism only (no divergence).
                        permut_pair = list(incodonssetlist)
                        basediff, inP, inD = CountBaseDiff(permut_pair), 1, 0
                        if basediff == 1 :
                            filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, permut_pair, inP, inD, filePs, filePn, fileDs, fileDn)
                        elif basediff > 1 :
                            int_permut_pairs = GetPairsInPermut(GetBestIntermediaryCodons(permut_pair, STOPcodons))
                            if int_permut_pairs != None :
                                for pair in int_permut_pairs :
                                    filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, pair, inP, inD, filePs, filePn, fileDs, fileDn)
                            else :
                                PrintWarningMessages(file, codon, "NP")
                        else :
                            PrintWarningMessages(file, codon, "UE")
                    elif len(set(codon_in_dict.values())) == 1 and len(set(codon_out_dict.values())) == 1 : #Divergence only (no polymorphism).
                        permut_pair = list(set(codon_in_dict.values()+codon_out_dict.values()))
                        basediff, inP, inD = CountBaseDiff(permut_pair), 0, 1
                        if basediff == 1 :
                            filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, permut_pair, inP, inD, filePs, filePn, fileDs, fileDn)
                        elif basediff > 1 :
                            int_permut_pairs = GetPairsInPermut(GetBestIntermediaryCodons(permut_pair, STOPcodons))
                            if int_permut_pairs != None :
                                for pair in int_permut_pairs :
                                    filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, pair, inP, inD, filePs, filePn, fileDs, fileDn)
                            else :
                                PrintWarningMessages(file, codon, "NP")
                        else :
                            PrintWarningMessages(file, codon, "UE")
                    else :
                        PrintWarningMessages(file, codon, "UE")
                elif len(list(set(codon_in_dict.values()+codon_out_dict.values()))) > 2 : #Checking for polymorphisms and divergences: more complex case when there is at least 3 possible codons in all sequences.
                    totcodons = list(set(codon_in_dict.values()+codon_out_dict.values()))
                    nbtotcodons = len(totcodons)
                    STOPdetected = STOPCount(STOPcodons, totcodonslist)
                    fileSTOP += STOPdetected
                    if len(list(set(codon_in_dict.values()+codon_out_dict.values())))+1 == len(incodonssetlist)+len(outcodonssetlist) : #Polymorphisms only (no divergence).
                        chosen_permut, outcodonatend = GetMostParcimoniousPermutation(GetScoredPermutations(totcodons), codon_out_dict)
                        permut_pairs = GetPairsInPermut(chosen_permut)
                        for permut_pair in permut_pairs :
                            basediff, inP, inD = CountBaseDiff(permut_pair), 1, 0
                            if basediff == 1 :
                                filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, permut_pair, inP, inD, filePs, filePn, fileDs, fileDn)
                            elif basediff > 1 :
                                int_permut_pairs = GetPairsInPermut(GetBestIntermediaryCodons(permut_pair, STOPcodons))
                                if int_permut_pairs != None :
                                    for pair in int_permut_pairs :
                                        filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, pair, inP, inD, filePs, filePn, fileDs, fileDn)
                                else :
                                    PrintWarningMessages(file, codon, "NP")
                    elif len(list(set(codon_in_dict.values()+codon_out_dict.values()))) == len(incodonssetlist)+len(outcodonssetlist) : #Polymorphisms with a divergence.
                        chosen_permut, outcodonatend = GetMostParcimoniousPermutation(GetScoredPermutations(totcodons), codon_out_dict)
                        permut_pairs = GetPairsInPermut(chosen_permut)
                        for permut_pair in permut_pairs :
                            basediff = CountBaseDiff(permut_pair)
                            inP, inD = 0, 0
                            if outcodonatend == 0 :
                                inP = 1
                            elif permut_pair[0] in incodonssetlist and permut_pair[-1] in incodonssetlist :
                                inP = 1
                            elif (permut_pair[0] in incodonssetlist and permut_pair[-1] in outcodonssetlist) or (permut_pair[-1] in incodonssetlist and permut_pair[0] in outcodonssetlist):
                                inD = 1
                            else :
                                PrintWarningMessages(file, codon, "UE")
                            if basediff == 1 :
                                filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, permut_pair, inP, inD, filePs, filePn, fileDs, fileDn)
                            elif basediff > 1 :
                                int_permut_pairs = GetPairsInPermut(GetBestIntermediaryCodons(permut_pair, STOPcodons))
                                if int_permut_pairs != None :
                                    for pair in int_permut_pairs :
                                        filePs, filePn, fileDs, fileDn = IncrementPsPnDsDn(gencode_dict, pair, inP, inD, filePs, filePn, fileDs, fileDn)
                                else :
                                    PrintWarningMessages(file, codon, "NP")
                    else :
                        PrintWarningMessages(file, codon, "UE")
                else :
                    PrintWarningMessages(file, codon, "UE")
            calculatedNI = CalculateNeutralityIndex(filePs, filePn, fileDs, fileDn)
            calculatedPnPs = CalculatePnPs(filePs, filePn)
            outputlist.append("\t".join([str(x) for x in [file[:14], filePs, filePn, fileDs, fileDn, fileSTOP, calculatedNI, calculatedPnPs, totSs, totNSs, lenseq, nboutgroups, nbingroups, "\n"]]))
        else :
            PrintWarningMessages(file, codon, "MD")
    # for i in range(1, len(outputlist)) :
        # outputlist[i] += "\t".join([str(x) for x in ["", totSs, totNSs, nboutgroups, nbingroups, "\n"]])
    outputfile.writelines(outputlist)
