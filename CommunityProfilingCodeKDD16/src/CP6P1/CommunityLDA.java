package CP6P1;
/*
 * 2015-12-14 ** Clean version parallel, divide dataset by words
 * For community profiling, sigmoid link function for user-user and tweet-tweet, polya-gamma data agumented & gibbs sampling
 */
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import common.FileUtil;
import common.Stopwords;

public class CommunityLDA {


	public static void main(String args[]) throws IOException {

		String dataset = /*"toy2"*/"CP";
		String foldername = "TMinput/";
		String indir = "G:/mmelab/Reserved/Users/lainey/adsc/data/"+dataset+"/"+foldername;
		String outdir = "data/"+dataset;
		String modelParas = indir+"/modelParametersmtc50.txt";	//parameter setting file
		String filelist = indir+"userlist.txt";				//user list file
		String retweetpath = indir+"retweet.txt";		//no retweet list		
		String dataDir = indir+"tweetsrfrdr/";				//tweets file
		String tempOutDir = outdir + "/TMoutput/";			//output folder
		String tRefIxpath = indir +"tRefIx.txt";		//tRefIx list
		String uRefIxpath = indir +"uRefIx.txt";		//tRefIx list
		String ugroupspath = indir +"ugroups.txt";
		
		(new File(tempOutDir)).mkdirs();

		//read parameters setting from file
		ArrayList<String> modelSettings = new ArrayList<String>();
		getModelPara(modelParas, modelSettings);
		int T = Integer.parseInt(modelSettings.get(0));
		int A = Integer.parseInt(modelSettings.get(1));
		float beta = Float.parseFloat(modelSettings.get(2));
		float rho = Float.parseFloat(modelSettings.get(3));
		int iteration = Integer.parseInt(modelSettings.get(4));
		int saveStep = Integer.parseInt(modelSettings.get(5));
		int saveTimes = Integer.parseInt(modelSettings.get(6));	
		T = A = 150;
		System.err.println("Topics:" + T + ", beta:" + beta + ", rho:" + rho + ", iteration:" + iteration + ", saveStep:" + saveStep
				+ ", saveTimes:" + saveTimes);
		modelSettings.clear();
		
		tempOutDir += ("C"+A+ "/");
		(new File(tempOutDir)).mkdirs();

		//read user list from file
		ArrayList<String> files = new ArrayList<String>();
		FileUtil.readLines(filelist, files);
		
		//read user list from file
//		ArrayList<String> retweet = new ArrayList<String>();
//		FileUtil.readLines(retweetpath, retweet);
		HashSet<String> retweet = new HashSet<String>();
		FileUtil.readLines(retweetpath, retweet);
		
		ArrayList<Integer> ltRefIx = new ArrayList<Integer>();
		FileUtil.readLines1(tRefIxpath, ltRefIx);
		
		ArrayList<Integer> luRefIx = new ArrayList<Integer>();
		FileUtil.readLines1(uRefIxpath, luRefIx);
		
		ArrayList<ArrayList<Integer>> ugroups = new ArrayList<ArrayList<Integer>>();
//		FileUtil.readLines2(ugroupspath, ugroups);
		ArrayList<Integer> lu = new ArrayList<Integer>();
		
		// the array to store all users
		ArrayList<User> docs = new ArrayList<User>();
		//user-user graph
		HashMap<Integer, List<Integer>> userNeighbour = new HashMap<Integer, List<Integer>>(); // user, the list of users who follow her, store the index of user, not the real user id, the index is the row number of each user in the file userlist.txt
		//user individual feature
		double[][] indfea = new double[files.size()][2];
		
		// 1. read users data from file
		new Stopwords();
		for (int i = 0; i < files.size(); i++) {
			System.out.println("Reading user "+i);
			User doc = new User(dataDir + files.get(i), i, files);
			docs.add(doc);
			lu.add(i);
		}
		ugroups.add(lu);
		// output wordMap and itemMap
		int wordMapSize = 2316020;
		int tsMapSize = 1633;
//		int wordMapSize = 9;
//		int tsMapSize = 4;
		
		//read graph
		Graphs g = new Graphs();
		String UserNetworkFileName = indir+"uunetwork.txt";
		userNeighbour = FileUtil.readHashMap(UserNetworkFileName);
		indfea = g.readIndFea(indir+"indfeature.txt", files.size(), 2);
//		indfea = g.readIndFea(indir+"indf2.txt", files.size(), 2);
		int para = 1, paranu = 1, parapopu = 1, parag=10;
		String outputDir = tempOutDir +"pg_" + parag +"_beta_"+beta+ "_rho_0.1/";
		FileUtil.mkdir(new File(outputDir));
//		FileUtil.mkdir(new File(outputDir + "/TaggedDocs/"));
		FileUtil.mkdir(new File(outputDir + "/modelRes/"));

		Model model = new Model(T, A, beta, rho,parag/*,para, paranu, parapopu*/);
		int initer = 0;  
		model.init(docs, wordMapSize, tsMapSize, userNeighbour, indfea,retweet,ltRefIx,luRefIx,ugroups);
		//model.init_frIter(docs, wordMapSize, tsMapSize,userNeighbour, indfea, outputDir+"modelRes/iter"+initer+"/model.",retweet,ltRefIx,luRefIx,ugroups);
//		model.printnoretweet(docs, outputDir + "/modelRes/");
		model.inference(iteration, docs, saveStep, saveTimes, outputDir + "/modelRes/",initer);
//		model.getResFromLastIteration(docs);
//		model.computeModelParameter();
//
//		// 3. output model results
//		// get uniWordMap and uniItemMap
//		System.out.println("saving the model...");
//		try {
//			model.saveModel(outputDir/*, uniWordMap*/);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		System.out.println("Ouputting tagged Docs...");
//		model.outTaggedDoc(docs, /*uniWordMap,*/ outputDir + "/TaggedDocs/");

		docs.clear();
		System.out.println("done");
		
	}

	public static void printDocs(ArrayList<User> docs) {
		for (int i = 0; i < docs.size(); i++) {
			for (int j = 0; j < docs.get(i).getDocWords().length; j++) {
				for (int m = 0; m < docs.get(i).getDocWords()[j].length; m++) {
					System.out.print(docs.get(i).getDocWords()[j][m] + "\t");
				}
				System.out.println("TS" + docs.get(i).getDocTimeStamp()[j] + "\t");
			}
		}
	}

	private static void getModelPara(String modelParas, ArrayList<String> modelSettings) {
		modelSettings.clear();
		// T , epsilon , beta , iteration , saveStep, saveTimes
		modelSettings.clear();
		modelSettings.add("100");
		modelSettings.add("100");
		modelSettings.add("0.01");
		modelSettings.add("1");
		modelSettings.add("1000");
		modelSettings.add("1000");
		modelSettings.add("1");

		ArrayList<String> inputlines = new ArrayList<String>();
		FileUtil.readLines(modelParas, inputlines);
		for (int i = 0; i < inputlines.size(); i++) {
			int index = inputlines.get(i).indexOf(":");
			String para = inputlines.get(i).substring(0, index).trim()
					.toLowerCase();
			String value = inputlines.get(i)
					.substring(index + 1, inputlines.get(i).length()).trim()
					.toLowerCase();
			switch (ModelParas.valueOf(para)) {
			case topics:
				modelSettings.set(0, value);
				break;
			case communities:
				modelSettings.set(1, value);
				break;
			case beta:
				modelSettings.set(2, value);
				break;
			case rho:
				modelSettings.set(3, value);
				break;
			case iteration:
				modelSettings.set(4, value);
				break;
			case savestep:
				modelSettings.set(5, value);
				break;
			case savetimes:
				modelSettings.set(6, value);
				break;
			default:
				break;
			}
		}
	}
	
	public enum ModelParas {
		topics,communities, beta, rho, iteration, savestep, savetimes;
	}

}
