package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;
import me.songbx.service.IfIntron;


public class CompareCNSldWithNoneCNSld {
	public static void main(String[] args) {
		for( int chrI=1; chrI<=10; chrI=chrI+1 ) {	
			try {
				
				IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
				
				PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSLDpattern/detail");
				
				String bwFilePosition =   "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw";
	    		Path bwFile = Paths.get(bwFilePosition);
	    		WigFileReader wig = WigFileReader.autodetect(bwFile);
	    		
	    		HashSet<String> checkedSnps = new HashSet<String>();
	    		
				String chr =  Integer.toString(chrI);
				
				int cnsLdCount = 0;
				int interLdCount = 0;
				int nonCnsLdCount = 0;
				
				double cnsLdSum = 0.0;
				double interLdSum = 0.0;
				double nonCnsLdSum = 0.0;
				int lineIndex = 0;
				
				File file1 = new File("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/chr"+Integer.toString(chrI)+"_ld.ld");
				BufferedReader reader1 = new BufferedReader(new FileReader(file1));
	            String tempString = reader1.readLine();
	            HashSet<String> cnsSnps = new HashSet<String>();
	            HashSet<String> cdsSnps = new HashSet<String>();
				while ((tempString = reader1.readLine()) != null) {
					++lineIndex;
					tempString = tempString.trim();
					String[] currencies = tempString.split("\\s+");
					chr=currencies[0];
					if( ! checkedSnps.contains(currencies[1]) ) {
						int position = Integer.parseInt(currencies[1]);
						Contig result = wig.query(chr, position, position);
						double thisMean = result.mean();
						if (thisMean > 0){
							cnsSnps.add(currencies[1]);
						}
						
						if(ifIntron.getElement(chr, position) == 2) {
							cdsSnps.add(currencies[1]);
						}
					}
					if( ! checkedSnps.contains(currencies[4]) ) {
						int position = Integer.parseInt(currencies[4]);
						Contig result = wig.query(chr, position, position);
						double thisMean = result.mean();
						if (thisMean > 0){
							cnsSnps.add(currencies[4]);
						}
						if(ifIntron.getElement(chr, position) == 2) {
							cdsSnps.add(currencies[4]);
						}
					}
					checkedSnps.add(currencies[1]);
					checkedSnps.add(currencies[4]);
					
					
					int numberOfCns = 0;
					if( cnsSnps.contains(currencies[1]) ) {
						numberOfCns = numberOfCns + 1;
					}
					if( cnsSnps.contains(currencies[4]) ) {
						numberOfCns = numberOfCns + 1;
					}
					
					int numberOfCds = 0;
					if( cdsSnps.contains(currencies[1]) ) {
						numberOfCds = numberOfCds + 1;
					}
					if( cdsSnps.contains(currencies[4]) ) {
						numberOfCds = numberOfCds + 1;
					}
					
					if( numberOfCns == 0 ) {
						nonCnsLdCount = nonCnsLdCount + 1;
						nonCnsLdSum = nonCnsLdSum + Double.parseDouble(currencies[6]);
						if( numberOfCds == 0 ) {
							outPut.write(currencies[6] + "\tnon-CNS\n");
						}else if( numberOfCds == 1 ) {
							outPut.write(currencies[6] + "\tCDS-boundary\n");
						}else {
							outPut.write(currencies[6] + "\tCDS\n");
						}
					}else if( numberOfCns == 1 ) {
						
						interLdCount = interLdCount + 1;
						interLdSum = interLdSum + Double.parseDouble(currencies[6]);
						
						if( numberOfCds == 0 ) {
							outPut.write(currencies[6] + "\tCNS-boundary\n");
						}else if( numberOfCds == 1 ) {
							outPut.write(currencies[6] + "\tCNS-CDS\n");
						}else {
							outPut.write(currencies[6] + "\tCNS-CDS-CDS\n");
						}
						
					}else if( numberOfCns == 2 ) {
						cnsLdCount = cnsLdCount + 1;
						cnsLdSum = cnsLdSum + Double.parseDouble(currencies[6]);
						if( numberOfCds == 0 ) {
							outPut.write(currencies[6] + "\tCNS\n");
						}else if( numberOfCds == 1 ) {
							outPut.write(currencies[6] + "\tCNS-CNS-CDS\n");
						}else {
							outPut.write(currencies[6] + "\tCNS-CNS-CDS-CDS\n");
						}
					}
					if( 0 == lineIndex%50000000) {
						System.out.println("lineIndex:"+lineIndex);
					}
				}
				reader1.close();
				
				System.out.println("cnsLdCount:" + cnsLdCount + " interLdCount:"+interLdCount+" nonCnsLdCount:"+nonCnsLdCount);
				System.out.println("cnsLdSum:" + cnsLdSum + " interLdSum:"+interLdSum+" nonCnsLdSum:"+nonCnsLdSum);
				
				outPut.close();
			} catch (IOException e) {
				e.printStackTrace();
			} catch (WigFileException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
    }
}
