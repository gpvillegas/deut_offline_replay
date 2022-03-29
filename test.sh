#!/bin/bash


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo
   echo "Syntax: scriptTemplate [-g|h|v|V]"
   echo "options:"
   echo "g     Print the GPL license notification."
   echo "h     Print this Help."
   echo "v     Verbose mode."
   echo "V     Print software version and exit."
   echo
}

############################################################

print_lic()
{
    echo "Printing GPL . . ."
}


############################################################
# Main program                                             #
############################################################
############################################################
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":hf:" option; do
    case $option in
	h) # display Help
            Help	    
	    exit;;       	
	f) # Enter a name
            fsys=$OPTARG;;	      
	\?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done


echo "Setting up $fsys filesystem directory on $HOSTNAME host"
