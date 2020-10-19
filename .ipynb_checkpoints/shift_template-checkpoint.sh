# Call as ./myScript -a -b /path/to/a/folder 

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
shift $((OPTIND - 1))  # shifts away every option argument,
                       # leaving your path as $1, and every
                       # positional argument as $@
path_argument="$1"
echo "$variable was chosen, path argument was $path_argument"