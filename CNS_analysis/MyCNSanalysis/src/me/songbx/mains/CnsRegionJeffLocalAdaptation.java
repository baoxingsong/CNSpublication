package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.service.IfIntron;

public class CnsRegionJeffLocalAdaptation {
	public static void readWigFile(String wigFile, HashMap<String, HashMap<Integer, Integer>> wigRecords) {
		try {
			File file = new File( wigFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String chr = "";
			Pattern p =Pattern.compile("chrom=(\\S+)");
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) == 't' ) {
					
				}else if( tempString.charAt(0) == 'v' ) {
					Matcher m = p.matcher(tempString);
					if( m.find() ) {
						chr = m.group(1);
					}
				}else {
					String[] arrOfStr = tempString.split("\\s+");
					if( wigRecords.containsKey(chr) ) {
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							if( wigRecords.get(chr).containsKey(Integer.parseInt(arrOfStr[0])) ) {
								wigRecords.get(chr).put(Integer.parseInt(arrOfStr[0]), wigRecords.get(chr).get(Integer.parseInt(arrOfStr[0]))+1 );
							}else {
								wigRecords.get(chr).put( Integer.parseInt(arrOfStr[0]), 1 );
							}
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		try {
			HashMap<String, HashMap<Integer, Integer>> cnsRecords = new HashMap<String, HashMap<Integer, Integer>> ();  // chr and positions
			HashMap<String, HashMap<Integer, Integer>> commonCnsVariantRegionCnsRecords = new HashMap<String, HashMap<Integer, Integer>> ();  // chr and positions
			HashMap<String, HashMap<Integer, Integer>> sorghumCnsRecords = new HashMap<String, HashMap<Integer, Integer>> ();  // chr and positions
			{
				for( int i=1; i<=10; ++i ) {
					String chr = Integer.toString(i);
					cnsRecords.put(chr, new HashMap<Integer, Integer>());
					commonCnsVariantRegionCnsRecords.put(chr, new HashMap<Integer, Integer>());
					sorghumCnsRecords.put(chr, new HashMap<Integer, Integer>());
				}
				System.out.println("wig file reading begin");
				readWigFile("/media/bs674/1_8t/AndCns/CNSBasedGwas/common_cns_variants.wig", commonCnsVariantRegionCnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5.wig", cnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5.wig", cnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.wig", cnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5.wig", cnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5.wig", cnsRecords);
				readWigFile("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.wig", sorghumCnsRecords);
			}
			
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/PCADAPT_pvals_depth.csv");
			for( int chro=1; chro<=10; ++chro ) {
				String chr = Integer.toString(chro);
				
				IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
				
				// here and allGoodPositions are callable and non-masked regions
				File file = new File("/media/bs674/1_8t/AndCns/overLapWithDomestication/localAdaptation/PCADAPT_pvals.csv");
				BufferedReader reader = new BufferedReader(new FileReader(file));
				String tempString = reader.readLine(); // skip the first line
				while ((tempString = reader.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+"); 
					if ( arrOfStr[0].compareTo(chr)==0 ){
						int position = Integer.parseInt(arrOfStr[1]);
						if ( commonCnsVariantRegionCnsRecords.get(chr).containsKey(position) ) {
							outPut.println(tempString+"\t-1\t"+"common-cns"+"\t" + ifIntron.getElement(chr, position)); // variants located in the common CNS variant region
						}else if( cnsRecords.get(chr).containsKey(position) ) {
							outPut.println(tempString+"\t" + cnsRecords.get(chr).get(position) +"\t"+"cns"+"\t" + ifIntron.getElement(chr, position));
						}else {
							outPut.println(tempString+"\t"+"0\tnon-cns"+"\t" + ifIntron.getElement(chr, position));
						}
						if ( sorghumCnsRecords.get(chr).containsKey(position) ) {
							outPut.println(tempString+"\t-2\t"+"sorghum-cns"+"\t" + ifIntron.getElement(chr, position)); // variants located in the common CNS variant region
						}
					}
				}
				reader.close();
			}
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
