cat hg38_to_GRCh37.chain | perl -lane '
    if($F[0] eq "chain"){
          if($F[2] ne $F[7] || $F[2] ne $ENV{"CHR"}){
                 $skip=1;
          }
          else{$curchr=$F[2];our $current=$F[5];our $coffset=$F[10]-$F[5];$skip=0;}
    }
elsif(!$skip) {
	my $offset=$F[2]-$F[1];
	print "$curchr $current ", $current+$F[0], " $coffset"; 
	$current=$current+$F[0]+$F[1];
	$coffset=$coffset+$offset;
}'
