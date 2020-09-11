package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;

public class Extractbambyreadsnames {
	public static void main(String[] args) {
		System.out.println("bam_file read_list output");
		
		try {
			PrintWriter outPut = new PrintWriter(args[2]);
			
			// put all the reads name in a hashset
			HashSet<String> readNames = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new FileReader(new File(args[1])));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				readNames.add(tempString.trim());
            }
			reader.close();
			
			// reads sam file one record by one record
			BufferedReader reader0 = new BufferedReader(new FileReader(new File(args[0])));
			tempString = null;
			while ((tempString = reader0.readLine()) != null) {
				tempString = tempString.trim();
				String[] arrOfStr = tempString.split("\\s+");
				if( readNames.contains(arrOfStr[0]) ) {
					outPut.println(tempString);
				}
            }
			reader.close();
			outPut.close();
		}catch (IOException e) {
            e.printStackTrace();
        }
	}
}
