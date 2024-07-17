import os

#----------------------------------------------------------
# F U N C T I O N   F I L E E X I S T S
#----------------------------------------------------------
# Determine if a file exists
# result = function(file)
#----------------------------------------------------------
def fileexists(file):
    if not os.path.exists(file):
        print(file + ' NOT FOUND.\tPROGRAM EXIT.')
        exit()
#----------------------------------------------------------
