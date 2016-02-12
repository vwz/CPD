package CP6P1;
import java.util.ArrayList;
import java.util.HashMap;
import common.FileUtil;


/**
 * A Document object represents a user, the document has been preprocessed, all words format as wordid and time format as timestampid
 */

public class User {
	
	// the ID or title of this document
	private int id;

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	private int [][] docWords;	//tweets 
	//use hashmap instead of two-d array for quick search when sample z, we need to find the retweet info for each tweet, need to be efficient
	HashMap<Integer, int[]> retweets;
	
	private int [] docTimestamp;
	
	// Init. Document -> get DataLines from Document
	public User(String Dir, int id, ArrayList<String> useridlist) {

		this.id = id;

		ArrayList<String> datalines = new ArrayList<String>();
		FileUtil.readLines(Dir, datalines);
		
		int totallinenb = datalines.size();

		this.docWords = new int[totallinenb][];
		this.docTimestamp = new int[totallinenb];
		ArrayList<int[]> tempretweet = new ArrayList<int[]>();		

		int lineNo = 0;
		for(int j=0; j<datalines.size(); j++) {
			String line = datalines.get(j);			
			String[] terms = line.split("\t");				
			String[] tweetbody = terms[2].split(" ");
			docWords[lineNo] = new int[tweetbody.length];
			for(int w = 0; w < tweetbody.length; w++)
				docWords[lineNo][w] = Integer.parseInt(tweetbody[w]);
			docTimestamp[lineNo] = Integer.parseInt(terms[1]);	
			if(terms[0].equals("R"))
			{
				int[] rpair = new int[]{lineNo,useridlist.indexOf(terms[terms.length-2]),Integer.parseInt(terms[terms.length-1])};
				tempretweet.add(rpair);
			}	
			lineNo++;
		}
		
		if(tempretweet.size()>0)
		{
			retweets = new HashMap<Integer, int[]>();
			for(int i=0; i<tempretweet.size();i++)
			{
				int[] pair = new int[]{tempretweet.get(i)[1],tempretweet.get(i)[2]}; 
				retweets.put(tempretweet.get(i)[0], pair);
			}
		}
		tempretweet.clear();
		datalines.clear();
	}
	

	public int[][] getDocWords() {
		return docWords;
	}

	public void setDocWords(int[][] docWords) {
		this.docWords = docWords;
	}

	public int[] getDocTimeStamp() {
		return docTimestamp;
	}

	public void setDocTimeStamp(int[] docts) {
		this.docTimestamp = docts;
	}

	public HashMap<Integer, int[]> getRetweets() {
		return retweets;
	}

	public void setRetweets(HashMap<Integer, int[]> retweets) {
		this.retweets = retweets;
	}

	
}