package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.model.BedRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.tribble.readers.TabixReader;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;

import java.nio.file.Path;
import java.nio.file.Paths;
/*
 * not done yet
 * */

public class CNSsnpDensityVsNonCNSRegion {
	public static void main(String[] args) {
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/overlapCNSwitheQTLAndSnp/panandCNS/snpdensity_java");
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46_cds.fa").getChromoSomeHashMap();
			
			
			HashMap<String, TabixReader> tabixReaders = new HashMap<String, TabixReader>();
			
			outPut.print("GenetypeID\tlength");
	        try {
	        	File file = new File("/media/bs674/and-CNS/maize282Genotyping/callablePipeline/callable_list");
	        	BufferedReader reader = new BufferedReader(new FileReader(file));
	            String tempString = null;
	            while ((tempString = reader.readLine()) != null) {
	            	tempString = tempString.trim();
	        		String baseName = tempString;
	        		baseName = baseName.replaceAll(".*\\/", "");
	        		baseName = baseName.replaceAll("_callable_status.bed.gz", "");
	        		outPut.print("\t" + baseName);
            		
					TabixReader tr = new TabixReader(tempString);
					tabixReaders.put(baseName, tr);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
	        
	        String bwFilePosition = "/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and_cns_depth.bw";
	        Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig0 = WigFileReader.autodetect(bwFile);
			
	        SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/panCnsAndCoreCns/pan_and.bam"));
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
//	            		Contig result = wig0.query(chr, i, i);
//	        			double thisMean = result.mean();
//	        			if( thisMean>0 ) {
//	        				totalLength = totalLength + 1;
//	        				nonMaskingPositions.add(i);
//	        			}
	            	}
	            }  
			}
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
	static boolean ifCallAble( int position, ArrayList<BedRecord> bedRecords, WigFileReader wig, String chr, int position2, WigFileReader wig2, String seq_name ) {
		for( BedRecord bedRecord : bedRecords ) {
			if( position >= bedRecord.getStart() && position <= bedRecord.getEnd() ) {
				return true;
			}
			if( bedRecord.getStart() > position ) {
				break;
			}
		}	
		Contig result;
		try {
			result = wig.query(chr, position, position);
			double thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
			
			result = wig2.query(seq_name, position2, position2);
			thisMean = result.mean();
			if (thisMean > 0){
				return true;
			}
		} catch (WigFileException | IOException e) {
			e.printStackTrace();
		}		
		return false;
	}
}

