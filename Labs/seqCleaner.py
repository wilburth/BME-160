# C:\Users\wwsch\anaconda3\python.exe
# Name: William Schlough(wschloug)
# Group Members: N/A

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example:
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAString (str):
    def length (self):
        return len (self) #changed length to len


    ''' 
    Return an uppercase version of the string, collapsing a single run of Ns.    
    '''
    def purify(self):
        upperData = self.upper() #make user input uppercase

        countN = str(upperData.count("N")) #count number of occurrences of "N" --> also use str() to format the number as a string value (joinData wont work unless every argument is of same type)

        positionN = upperData.find("N") #find first occurrence of "N"

        startData = upperData[ : positionN] #from start up to but not including positionN
        endData = upperData[positionN + 1 : ] #from positionN + 1 (omitting positionN)

        joinData = startData + "{" + countN + "}" + endData #join the two sections of the original string (minus the first occurrence of "N"), along with the original number of occurrences of "N"
        cleanData = joinData.replace("N", "") #replace all remaining "N" with "" (effectively removing them)

        return cleanData #this returns the final version of the string, formatted properly (this is a string obj with data "startData + { + countN + } + endData")


def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?') #any input taken is as a string
    thisDNA = DNAString (data) #initializes data as class DNAString (makes it an str value)
    pureData = thisDNA.purify() #uses function purify on data within thisDNA --> pureData is filled with this
    print (pureData)

main()
