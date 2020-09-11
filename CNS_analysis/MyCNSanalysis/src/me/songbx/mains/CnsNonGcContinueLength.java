package me.songbx.mains;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.GetContinueNonGcs;
import me.songbx.model.ContinueNonGc;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CnsNonGcContinueLength {
	public static void main(String[] args) {
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		
		new CnsNonGcContinueLength(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsIntronUtrCns.sam", "/media/bs674/1_8t/AndCns/nongcDistribution/nonOpTfbsIntronUtrCns");
		new CnsNonGcContinueLength(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", "/media/bs674/1_8t/AndCns/nongcDistribution/nonIntronUtr_OpOrTFBS");
		
		new CnsNonGcContinueLength(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_highMethylation.bam", "/media/bs674/1_8t/AndCns/nongcDistribution/nonIntronUtrCns_highMethylation");
		new CnsNonGcContinueLength(maizeGenome, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_lowMethylation.bam", "/media/bs674/1_8t/AndCns/nongcDistribution/nonIntronUtrCns_lowMethylation");
	}
	public CnsNonGcContinueLength( HashMap<String, String> maizeGenome, String bamFile, String outputFile ) {
		int extendDistand = 200;
		try {
			HashMap<String, String> cnsSeqs =  new HashMap<String, String> ();
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart() - extendDistand;
	            if( start < 1 ) {
	            	start = 1;
	            }
	            int end = samRecord.getAlignmentEnd() + extendDistand;
	            String chr = samRecord.getContig();
	            if( end > maizeGenome.get(chr).length() ) {
	            	end = maizeGenome.get(chr).length();
	            }
	            String seq_name = chr+":"+start+"-"+end;
	            String seq = maizeGenome.get(chr).substring(start-1, end);
	            cnsSeqs.put(seq_name, seq);
			}
			reader0.close();
			PrintWriter outPut = new PrintWriter(outputFile);
			for( String seq_name : cnsSeqs.keySet()) {
				String seq = cnsSeqs.get(seq_name);
				ArrayList<ContinueNonGc>  continueNonGcs = GetContinueNonGcs.getContinueNonGcs( seq );
				for ( ContinueNonGc continueNonGc : continueNonGcs ) {
					outPut.println( seq_name + "\t" + seq.length() + "\t" + continueNonGc.getPosition() + "\t" + continueNonGc.getLength() );
				}
			}
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
