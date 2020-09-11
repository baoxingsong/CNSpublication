package me.songbx.mains;

import java.io.File;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.service.IfIntron;

public class OverLapCnsWithSorghumCds {
	public static void main(String[] args) {
		String bamFile;
		String outputFile;
		String gffFile;
		/*
		bamFile  = "/media/bs674/1_8t/AndCns/CNSprofiles/sorghum_nonOpTfbsEqtlCns.sam";
		outputFile = "/media/bs674/1_8t/AndCns/CNSprofiles/sorghum_nonOpTfbsEqtlCns_intronInMaize_cdsInSorghum_gff.sam";
		gffFile = "/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3";
		new OverLapCnsWithSorghumCds(bamFile, outputFile, gffFile);
		
		bamFile  = "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam";
		outputFile = "/media/bs674/1_8t/AndCns/CNSprofiles/sorghum_intronInMaize_cdsInSorghum_gff.sam";
		new OverLapCnsWithSorghumCds(bamFile, outputFile, gffFile);
		*/
		bamFile  = "/media/bs674/1_8t/AndCns/CNSprofiles/Setaria_nonOpTfbsEqtlCns.sam";
		outputFile = "/media/bs674/1_8t/AndCns/CNSprofiles/Setaria_nonOpTfbsEqtlCns_intronInMaize_cdsInSetaria_gff.sam";
		gffFile = "/media/bs674/2t/genomeSequence/Setaria_italica/Setaria_italica.Setaria_italica_v2.0.42.gff3";
		new OverLapCnsWithSorghumCds(bamFile, outputFile, gffFile);
		/*
		bamFile  = "/media/bs674/1_8t/AndCns/Setaria_italica/result/5.bam";
		outputFile = "/media/bs674/1_8t/AndCns/CNSprofiles/Setaria_intronInMaize_cdsInSetaria_gff.sam";
		new OverLapCnsWithSorghumCds(bamFile, outputFile, gffFile);
		*/
    }
	
	public OverLapCnsWithSorghumCds(String bamFile, String outputFile, String gffFile) {
		IfIntron ifIntron = new IfIntron(gffFile);
		try {
			PrintWriter outPut = new PrintWriter(outputFile);
			Pattern sorghumStartPositionPattern = Pattern.compile("^(\\d+)H");
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            Matcher m = sorghumStartPositionPattern.matcher(samRecord.getCigarString());
	            if( m.find() ) {
	            	int sorghumStart = Integer.parseInt(m.group(1));
	            	String sorghumSequence = samRecord.getReadString();
	            	if( samRecord.getReadNegativeStrandFlag() ) {
	            		sorghumSequence = ChromoSomeReadImpl.getReversecomplementary(sorghumSequence);
	            	}
	            	int overLapedWithCdsLength = 0;
	            	for( int sequenceIndex=0; sequenceIndex < sorghumSequence.length(); ++sequenceIndex ) {
	            		if( sorghumSequence.charAt(sequenceIndex) != 'n' &&  sorghumSequence.charAt(sequenceIndex) != 'N' ) {
	            			if(ifIntron.getElement(samRecord.getReadName(), sorghumStart+sequenceIndex)==2 ){
	            				++overLapedWithCdsLength;
	            			}
	            			if(overLapedWithCdsLength>=10) {
		            			outPut.print( samRecord.getSAMString() );
		            			break;
		            		}
	            		}
	            	}
	            }
			}
			reader0.close();
			outPut.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
