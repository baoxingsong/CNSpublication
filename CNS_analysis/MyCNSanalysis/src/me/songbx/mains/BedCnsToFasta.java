package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import me.songbx.impl.ChromoSomeReadImpl;

public class BedCnsToFasta {
	public static void main(String[] args) {
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
		new BedCnsToFasta(maizeGenome, "/media/bs674/1_8t/AndCns/coreCNSBasedGwas/core_cns.bed", "/media/bs674/1_8t/AndCns/coreCNSBasedGwas/core_cns.fa");
	}
	
	public BedCnsToFasta( HashMap<String, String> maizeGenome, String bedfile, String outputFile ) {
		try {
			HashMap<String, String> cnsSeqs =  new HashMap<String, String> ();
			
			try {
	        	File file = new File(bedfile);
	    		BufferedReader reader = new BufferedReader(new FileReader(file));
	            String tempString = null;
				while ((tempString = reader.readLine()) != null) {
					String[] arrOfStr = (tempString.trim()).split("\\s+");
					int start = Integer.parseInt(arrOfStr[1]) + 1;
		            int end = Integer.parseInt(arrOfStr[2]);
		            if ( (end - start) > 6 ) {
			            String chr = arrOfStr[0];
			            String seq_name = chr+":"+start+"-"+end;
			            String seq = maizeGenome.get(chr).substring(start-1, end);
			            cnsSeqs.put(seq_name, seq);
		            }
				}
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
			
			int cnsIndex = 0;
			PrintWriter outPut = new PrintWriter(outputFile);
			for( String seq_name : cnsSeqs.keySet()) {
				String seq = cnsSeqs.get(seq_name);
	            outPut.println(">" + cnsIndex + "\t" + seq_name);
	            outPut.println(seq);
	            cnsIndex++;
			}
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}

