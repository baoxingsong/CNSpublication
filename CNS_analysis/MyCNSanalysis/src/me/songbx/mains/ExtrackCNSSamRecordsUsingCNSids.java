package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.impl.ChromoSomeReadImpl;

public class ExtrackCNSSamRecordsUsingCNSids {
	public static void main( String args[] ) {
		new ExtrackCNSSamRecordsUsingCNSids("/media/bs674/1_8t/AndCns/CNSBasedGwas/common_cns_variants", "/media/bs674/1_8t/AndCns/CNSBasedGwas/common_cns_variants.sam");
		// extract the common CNS variants (MAF> i.e. 0.1) and output to a sam fil
	}
	public ExtrackCNSSamRecordsUsingCNSids(String cnsIdFile, String outPutFile) {
		try {
			HashSet<String> commonCnsIds = new HashSet<String>();
			File file = new File( cnsIdFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				commonCnsIds.add(arrOfStr[0]);
			}
			reader.close();
			
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
			String bwFilePosition = "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw";
			Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig0 = WigFileReader.autodetect(bwFile);
    		
			PrintWriter outPut = new PrintWriter(outPutFile);
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam"));
	        SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            String chr = samRecord.getContig();
	            HashSet<Integer> nonMaskingPositions = new HashSet<Integer>();
	        
	            int totalLength = 0;
	            for( int i=start-1; i<end; ++i ) {
	            	if( maizeGenome.get(chr).charAt(i) != 'n' ) {
	            		Contig result = wig0.query(chr, i+1, i+1);
	        			double thisMean = result.mean();
	        			if( thisMean>0 ) {
	        				totalLength = totalLength + 1;
	        				nonMaskingPositions.add(i);
	        			}
	            	}
	            }
	            //if(totalLength > 30) {
		            String seq_name = chr+":"+start+"-"+end;
		            if( commonCnsIds.contains(seq_name) ) {
		            	outPut.print(samRecord.getSAMString());
		            }
	            //}
			}
			outPut.close();
			wig0.close();
		} catch (IOException | WigFileException e) {
			e.printStackTrace();
		}
	}
}
