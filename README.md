[![Build Status](https://travis-ci.org/jokergoo/circlize.svg)](https://travis-ci.org/jokergoo/circlize)


## NovoSV: Identify and parse the pattern of chromosomal structural variation

NovoSV is an efficient tool for chromosomal variation detection, which can accurately identify abnormal connections and parse the pattern of chromosomal variations. NovoSV is open source PERL program and has been validated on GNU/Linux systems.

### INSTALL
NovoSV is written in perl, users can run run this program directly after download source code

However, NovoSV need the following tools installed.

- `bwa (above v.0.7.8)`
- `cap3`
- `samtools (above v.0.1.19)`

NovoSV intend to compatible with Linux systems. Though perl code can be decode on windows platform, most of the tools called by NovoSV are not compatible with windows platform.

### Running NovoSV
Instance of running NovoSV can be found is "sample" directory.

NovoSV take bam file, reference file as input. Users also needs specific the path of samtools, bwa and CAP3, we encourage users add the path of these tools to environment variables to simplify the input.

Users can also specific parameters that influence the behaviours of NovoSV, including thresholds of PE reads or soft-clipped reads and the properly attributes of clusters. For details parameters, see "perl NovoSV.pl -h".

### Result format
NovoSV will report both identified abnormal connections and structure variations. 

####Abnormal connections(ACs) can be found in "sample_name"_out.txt or user specified file.
ACs starts with a '>' symbol, end at the start of next ACs or the end of this file. Example of AC are shown in following lines:

>\>2	8	21310256	21310256	15	36	0.20	8	21309659	21309659	r1[r2[	34	0   13	2

>ST-E00142:35:H04DDALXX:2:1211:32272:52854 161 8 21310258 60 118M32S 21310258 = 21309715 -278 TAAATGGGCACAAAAATATTTGAGGTAGCATTGTTTGTAGAGCAAAAATTACAAACAACCTGATATTCATCACCAAGAGAACAGAAATATAAATTGCATATTTATACAATGAAAACTTGGTGCAATCACATTGGAGACCAATCTGCCACT >=>>??><??>?????????????>???9????6>????????????????>@?@>@@?>A;@@;?@<@>-@>-9@<?8@@,@??@@>@7@@?@<@@=<?@?3@?@7=?<,@?@A,>?@4A@,@6<.:8A@@8*A9A<9>>@.B8=79;; SA:Z:8,21309657,+,116S34M,12,0; MC:Z:150M BD:Z:JJLDOONKPPIIKDDDKLLKCKNKMMNNPOLKKILCKINNKKPOKDDDKKLMHKDKHKKMNNNMLLKJMLMMHMLKLKKKJKHNKJDKLLLLDKKKOOLLLKCLLLMHKKMOKEELNLLLNJPPLLNNIIMLLMLMMNONMNPMPPNLHM MD:Z:118 PG:Z:MarkDuplicates RG:Z:EC1120_H_H04DDALXX_L2 BI:Z:LLMGOOMKNOJKMFFFLLNKFMNKMMNMPNMKMJKFMJNMKKPNLFFFLKMMJLFLKLLLMMNMLNLLOMOOJMMMNLLLLMKOLLGMMOMMGMLNNONMOLGNOMNJLMMOMHHNPMOMOLOPNNPPKLPNPNNNNOPPPPRNOONMJO NM:i:0 MQ:i:60 AS:i:118 XS:i:0

>......

>PE

>ST-E00142:35:H04DDALXX:4:1104:15849:25640 161 8 21310258 60 115M35S 21310258 = 21309657 -384 ATGGGCACAAAAATATTTGAGGTAGCATTGTTTGTAGAGCAAAAATTACAAACAACCTGATATTCATCACCAAGAGAACAGAAATATAAATTGCATATTTATACAATGAAAACTTGGTGCAATCACATTGGAGACCAATCTGCCACTATA >>?>?@@?@????@@@??@@@?@@@?@=?@><=@@@@?@@@>>??@?@?@?@@@?=@B<@AAA@@AA>=@@A@<?@>@?A?A@@AA8A@@A@A<AAAA:?>A?6A@AA?5@?;A<B?=@BBA=@AAB9A@AA=@<4=A@>AA=@B=?;:8 SA:Z:8,21309657,+,113S37M,30,0; MC:Z:45S105M BD:Z:JJLKMQRJKLEEDKLLJCLNKMMNNPOLJLILCLINNKKPOKDDDKJLMIKDKIKKMMNNMLLJJMLMMHMLKLKKKJKINKJDKLLLLDKJLOOLLLJCLLLMIKKLNJDDLNLMLNJPPLLNNIJMKMLKLMNONMMOMPRPNIMNLL MD:Z:115 PG:Z:MarkDuplicates RG:Z:EC1120_H_H04DDALXX_L4 BI:Z:LLMLLPQKKMGGGLLNKFMOKMMNMPOMKMJKFMJNMKKPOLFFFLKMNJLGLJLLMNMPNLOLKOMOOIMMLMLLLLMKOLLGMMOMMGMLNNPNMOLGNOMNKMMMPLHHNOLOMOLOQNNPPKLOMOMNNNOOOOOROQQPNKOOOM NM:i:0 MQ:i:60 AS:i:115 XS:i:0

>......

- `'\>2': id of this ACs, ids are  to distinguish different ACs;`
- `'8	21310256	21310256': specific the left breakpoint of this AC;`
- `'15': number of soft-clipped reads support this AC;`
- `'36': read-depth at left breakpoint;`
- `'0.2': Confidence of a true AC;`
- `'8	21309659	21309659': right breakpoint;`
- `'r1[r2[': connection type coding in VCF format, r1 and r2 stand for left breakpoint and right breakpoint separately;`
- `'34': PE support reads;`
- `'0': Number of PE reads that inconsistent with soft-clipped reads;`
- `'13': number of soft-clipped reads that are detected directly by alignment information;`
- `'2': number of soft-clipped reads that are detected indirectly by reads sequence similarity;`

- `Following lines are soft-clipped support reads and PE support reads. This lines are written in a bam-like format, actually most informations are identical with bam file, except the forth column which is the location of breakpoint rather than the alignment of this reads.;`
- `The number of soft-clipped reads and PE support reads in these lines are consistent with previously described informations;`

####Structure variations results can be found in "sample_name"_struc.txt or user specified file.
each results start with a '>' symbol, end at the start of next SV results or the end of the this file. Example of SV results are shown in the following lines:

>\>0	3	translocation_copy	8:21310256	8:21309659	8:21309875	216	+/+

>8:21310256	r1[r2[	8:21309659	15	36	0.20	34	
8:21309659	]r2]r1	8:21310256	13	73	0.00	35	
8:21309875	r1[r2[	8:21310260	8	47	0.00	47	

- `'>0': id of this SVs, ids are  to distinguish different SVs;`
- `'3': the number of ACs that support this SV;`
- `'translocation_copy': specific the type of this SV, other type include "translocation_cut", "deletion", "CNV";`
- `'8:21310256': the insert size of this SV;`
- `'8:21309659' and '8:21309875' specific the original of the insert sequence;`
- `'216': the length of insert sequences;`
- `'+/+': the connection type of this SV;`

Following lines are ACs that support this SV, each column means:
- `the left breakpoint;`
- `the connection type between left breakpoint and right breakpoint;`
- `the right breakpoint;`
- `number of soft-clipped support reads;`
- `reads depth at left breakpoint;`
- `Confidence of a true AC;`
- `Number of soft-clipped reads that are detected directly by alignment information;`



