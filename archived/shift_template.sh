# Call as ./myScript /path/to/a/folder -a -b

path_argument="$1"
shift   # Shifts away one argument by default

while getopts ":ab" opt; do
 case $opt in
  a) 
   variable=a
   ;;
  b)    
   variable=b
   ;;
  \?)
   echo "invalid option -$OPTARG"
   exit 0
 esac
done

echo "$variable was chosen, path argument was $path_argument"