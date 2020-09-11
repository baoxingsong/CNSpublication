package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;

import edu.unc.genomics.Contig;
import edu.unc.genomics.io.WigFileException;
import edu.unc.genomics.io.WigFileReader;
import me.songbx.service.IfIntron;

public class ListCNSandNonCNS {
	public static void main(String[] args) {
		try {
			IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
			
			PrintWriter outPut0 = new PrintWriter("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/CDS_SNPlist");
			PrintWriter outPut1 = new PrintWriter("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/CNS_SNPlist");
			PrintWriter outPut2 = new PrintWriter("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/other_SNPlist");
			
			String bwFilePosition =   "/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bw";
    		Path bwFile = Paths.get(bwFilePosition);
    		WigFileReader wig = WigFileReader.autodetect(bwFile);
    		
			File file1 = new File("/media/bs674/2019junehackatho/282V4/uplifted_APGv4/maizeagp4.bim");
			BufferedReader reader1 = new BufferedReader(new FileReader(file1));
            String tempString = reader1.readLine();
			while ((tempString = reader1.readLine()) != null) {
				tempString = tempString.trim();
				String[] currencies = tempString.split("\\s+");
				String chr = currencies[0];
				int position = Integer.parseInt(currencies[3]);
				
				if(ifIntron.getElement(chr, position) == 2) {
					outPut0.write(currencies[1] + "\n");
				}else {
					Contig result = wig.query(chr, position, position);
					double thisMean = result.mean();
					if (thisMean > 0){
						outPut1.write(currencies[1] + "\n");
					}else {
						outPut2.write(currencies[1] + "\n");
					}
				}
			}
			reader1.close();
			
			outPut0.close();
			outPut1.close();
			outPut2.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (WigFileException e) {
			e.printStackTrace();
		}
    }
}
