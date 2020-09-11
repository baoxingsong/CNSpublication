package me.songbx.mains;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import me.songbx.impl.ChromoSomeReadImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SamCnsToFasta {
	public static void main(String[] args) {
		
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/wd/cnsLength/masked_B73_v4_k20_46_cds.fa").getChromoSomeHashMap();
		
		//new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam", "/media/bs674/1_8t/AndCns/mapCnsToReads/sorghum/cns_seq.fa");
		/*
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsIntronUtrCns.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsIntronUtrCns.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.fa");
		
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsIntronUtrCns.sam", "/media/bs674/1_8t/AndCns/orthorfinderClusterCns/knownUnknown/nonOpTfbsIntronUtrCns.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", "/media/bs674/1_8t/AndCns/orthorfinderClusterCns/knownUnknown/nonIntronUtr_OpOrTFBS.fa");
		
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_highMethylation.bam", "/media/bs674/1_8t/AndCns/orthorfinderClusterCns/highMethyLationLowMethylation/nonIntronUtrCns_highMethylation.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_lowMethylation.bam", "/media/bs674/1_8t/AndCns/orthorfinderClusterCns/highMethyLationLowMethylation/nonIntronUtrCns_lowMethylation.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntron_openchromatin.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntron_openchromatin.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_noneopenchromatin.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_noneopenchromatin.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow.fa");
		*/
		
		/*
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/loop.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/loop.fa");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.fa");
		*/
		
		new SamCnsToFasta(maizeGenome, "/media/bs674/1twdblack/allCNSbamFiles/sorghum.bam", "/media/bs674/wd/cnsLength/sorghum_cns.fasta", "/media/bs674/wd/cnsLength/sorghum_cns.length");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1twdblack/allCNSbamFiles/sugarcane_tareploid.bam", "/media/bs674/wd/cnsLength/sugarcane_tareploid_cns.fasta", "/media/bs674/wd/cnsLength/sugarcane_tareploid_cns.length");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1twdblack/allCNSbamFiles/Hyparrhenia_diplandra.bam", "/media/bs674/wd/cnsLength/Hyparrhenia_diplandra_cns.fasta", "/media/bs674/wd/cnsLength/Hyparrhenia_diplandra_cns.length");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1twdblack/allCNSbamFiles/Miscanthus_sinensis.bam", "/media/bs674/wd/cnsLength/Miscanthus_sinensis_cns.fasta", "/media/bs674/wd/cnsLength/Miscanthus_sinensis_cns.length");
		new SamCnsToFasta(maizeGenome, "/media/bs674/1twdblack/allCNSbamFiles/Chrysopogonserrulatus.bam", "/media/bs674/wd/cnsLength/Chrysopogonserrulatus_cns.fasta", "/media/bs674/wd/cnsLength/Chrysopogonserrulatus_cns.length");
		
	}
	
	
	public SamCnsToFasta( HashMap<String, String> maizeGenome, String bamFile, String outputFile, String outputFileLength ) {
		try {
			HashMap<String, String> cnsSeqs =  new HashMap<String, String> ();
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            String chr = samRecord.getContig();
	            String seq_name = chr+":"+start+"-"+end;
	            String seq = maizeGenome.get(chr).substring(start-1, end);
	            cnsSeqs.put(seq_name, seq);
			}
			reader0.close();
			int cnsIndex = 0;
			PrintWriter outPut = new PrintWriter(outputFile);
			PrintWriter outPut2 = new PrintWriter(outputFileLength);
			for( String seq_name : cnsSeqs.keySet()) {
				String seq = cnsSeqs.get(seq_name);
	            outPut.println(">" + cnsIndex + "\t" + seq_name);
	            outPut.println(seq);
	            cnsIndex++;
	            
	            seq = seq.replaceAll("n", "");
	            seq = seq.replaceAll("N", "");
	            outPut2.println(">" + cnsIndex + "\t" + seq_name + "\t" + seq.length());
			}
			outPut.close();
			outPut2.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}


/*
grep -v ">" /media/bs674/1_8t/AndCns/orthorfinderClusterCns/knownUnknown/nonOpTfbsIntronUtrCns.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'

n=649134
a=45802701
c=34525215
g=34750110
t=4653204
sum(c,g)/sum(a,t,c,g)
0.5785903

grep -v ">" /media/bs674/1_8t/AndCns/orthorfinderClusterCns/knownUnknown/nonIntronUtr_OpOrTFBS.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n=7470
a=40101109
c=42076616
g=42166806
t=40119892
sum(c,g)/sum(a,t,c,g)
0.5122288

grep -v ">" /media/bs674/1_8t/AndCns/orthorfinderClusterCns/highMethyLationLowMethylation/nonIntronUtrCns_highMethylation.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n 372351
a=25758629
c=18563991
g=18721886
t=26011765
sum(c,g)/sum(a,t,c,g)
0.4186777

grep -v ">" /media/bs674/1_8t/AndCns/orthorfinderClusterCns/highMethyLationLowMethylation/nonIntronUtrCns_lowMethylation.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n=284253
a=60283470
c=58124063
t=60779513
g=58278838
sum(c,g)/sum(a,t,c,g)
0.4901879

grep -v ">" /media/bs674/1_8t/AndCns/CNSprofiles/nonIntron_openchromatin.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n=7470
a=10522290
c=13978325
t=10533965
g=14048009

sum(c,g)/sum(a,t,c,g)
#0.5710036

grep -v ">" /media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n=510886
a=28319019
c=22448447
g=22589896
t=28992046
sum(c,g)/sum(a,t,c,g)
0.440045



grep -v ">" unknownCNS.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n 531990
a=30200156
c=23818985
g=23975322
t=30906656
sum(c,g)/sum(a,t,c,g)
0.438878

grep -v ">" loop.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n 754
a=10023531
c=6487690
g=6531281
t=10065901
sum(c,g)/sum(a,t,c,g)
0.3932226

grep -v ">" nonIntronUtr_OpOrTFBS.fa | awk -vFS="" '{for(i=1;i<=NF;i++)w[tolower($i)]++}END{for(i in w) print i,w[i]}'
n 7470
a=40101109
c=42076616
g=42166806
t=40119892
sum(c,g)/sum(a,t,c,g)



**/
