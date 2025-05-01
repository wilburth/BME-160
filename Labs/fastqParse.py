# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A


'''
Read the seqname line of a FASTQ file and displays each field's information on a new line for easy readability.
example input:
    @EAS139:136:FC706VJ:2:2104:15343:197393
'''

class FastqString (str):
    ''' class FastqString .'''


    def parse(self):
        '''
        parse out each field of the input string and returns a string with each element on a new line
        '''



        separateFastq = self.split(':') # splits at every occurrence of ':' (effectively creating list of 7 elements --> (0, 1, 2, 3, 4, 5, 6))

        instrument = separateFastq[0] #use the first element of separateFastq
        positionAtSign = instrument.split('@') #split instrument at the '@' sign --> create list of 2 elements
        instrumentFinal = positionAtSign[1] #use the 2nd element in list positionAtSign (this will be the element that follows the '@')

        runId = separateFastq[1] #use the 2nd element in list separateFastq
        flowCellId = separateFastq[2] #use the 3rd element in list separateFastq
        flowCellLane = separateFastq[3] #use the 4th element in list separateFastq
        titleNumber = separateFastq[4] #use the 5th element in list separateFastq
        xCoord = separateFastq[5] #use the 6th element in list separateFastq
        yCoord = separateFastq[6] #use the 7th element in list separateFastq

        pureFastq = 'Instrument: ' + instrumentFinal + '\n' + 'Run ID: ' + runId + '\n' + 'Flow Cell ID: ' + flowCellId + '\n' + 'Flow Cell Line: ' + flowCellLane + '\n' + 'Title Number: ' + titleNumber + '\n' + 'X-coord: ' + xCoord + '\n' + 'Y-coord: ' + yCoord

        return pureFastq


def main():
    '''
    take an input string of 7 fields, sets it to be of class FastqString, calls method 'parse' on it, and prints it.
    '''

    fastqInput = input('FASTQ seqname: ') #any input taken is as a string
    thisFASTQ = FastqString (fastqInput) #initializes data as class FastqString
    pureFASTQ  = thisFASTQ.parse() #uses function parse on data within thisFASTQ --> pureFASTQ is filled with this
    print(pureFASTQ)

main()
