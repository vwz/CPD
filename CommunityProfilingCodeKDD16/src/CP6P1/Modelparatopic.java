package CP6P1;
/*
 * Community profiling - Hongyun Cai & 14/09/2015
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import Distribution.PolyaGamma;
import common.ComUtil;
import common.FileUtil;
import common.MatrixUtil;

public class Modelparatopic {
	// user related parameters
	// public Document Doc;
	public int T; // no of topics
	public int U; // no of users
	public int V; // vocabulary size
	public int M; // no of time stamps
	public int A; //no of communities
	public int Nf = 2;	//no of features for each user
	public int N2f = 2*Nf;
	public int NLu = 3589811;
	public int NLt =  992522;
//	public int NLu = 8;
//	public int NLt =  8;
	public double para = 2;
	public double paranu = 3;
	public double parapopu = 0.5;

	// hyperparameters
	public float alpha;
	public float beta;	
	public float rho;


	// model parameters   
	public int niters; // number of Gibbs sampling iteration

	// Estimated/Inferenced parameters
	public float[][] theta; // community-topic distribution, C*T
	public float[][] vPhi; // topic-word distribution, T*V
	public float[][] pai; //user-community distribution, U x C
	public double[][][] eta; //community-community retweet link distribution on topic k,  T x C x C 
	public double[] nu; // individual factor  
	//augmented variable
	public float[] lambda;	// U x U (any two users), NLu
	public float[] delta;	//N x N (any two tweets), NLt

	// Temp variables while sampling
	public int[][] Z; // U x N
	public int[][] C; //U x N
	public int[][] NTW; // T x V, sum is: SNTW[T]
	public int[][] NUC; // sum U x C, sum is: SNUC[U]; c_hat before normalisation
	public int[][] NCT; //C x T, sum is SNCT[C]
	public int[][] popularity;	//T x M, topic popularity at timeslice t
	
	public float[] SNTW;
	public float[] SNUC;
	public float[] SNCT;
	public float[] SNTC;	//for normalize hat_z_k
	public float[] SPTM;	//for normalize popularity
	
	//intermediate variables for training
	public PolyaGamma[] m_uPGsampler;
	public PolyaGamma[] m_tPGsampler;
	
	//variable record the neighbour of each user and tweet   ?????????
	HashMap<Integer, List<Integer>> userNeighbour; // user, the list of users who follow her, store the index of user, not the real user id, the index is the row number of each user in the file userlist.txt
	HashMap<Integer, List<int[]>> userInNeighbours;	//users and the list of users who she follows, also the index of userNeighbour to ref lambda
 	HashMap<String, List<int[]>> tweetInNeighbours; //the tweet and the list of tweets which retweet it, for each int[], we store the user and the index of the tweet in the user's tweets, also for the second int[], we stire the ubdex to refdelta
	HashSet<String> retweet;	
 	
	//individual statistics
	double[][] ind_fea;	//U * number of features; feature need normalize to keep in the same scale
	
	public double[] ptopic;
	
	public Modelparatopic(int no, int comno, float betaV, float rhoV) {
		T = no;
		A = comno;
		beta = betaV;
		//alpha = (float) (50.0 / T);
		alpha = 1f;
		rho = rhoV;
	}

	/**
	 * initialize the model
	 */
	protected boolean init(ArrayList<User> users, int vocabSize, int tsNo,HashMap<Integer, List<Integer>> myuserNeighbour, double[][] myindfea,HashSet<String> myretweet) {
		U = users.size();
		V = vocabSize;
		M = tsNo;
		// assign topics and communities randomly
		Z = new int[U][];
		C = new int[U][];		
		for (int i = 0; i < U; i++) {
			Z[i] = new int[users.get(i).getDocWords().length];
			C[i] = new int[Z[i].length];
			for (int j = 0; j < Z[i].length; j++) {
 				Z[i][j] = (int) Math.floor(Math.random() * T);
				if (Z[i][j] < 0)
					Z[i][j] = 0;
				if (Z[i][j] > T - 1)
					Z[i][j] = (int) (T - 1);
				C[i][j] = (int) Math.floor(Math.random() * A);
				if (C[i][j] < 0)
					C[i][j] = 0;
				if (C[i][j] > A - 1)
					C[i][j] = (int) (A - 1);				
			}
		}
		m_uPGsampler = new PolyaGamma[NLu];
		for(int i=0; i<NLu; i++)
		{
			m_uPGsampler[i] = new PolyaGamma();
		}
		m_tPGsampler = new PolyaGamma[NLt];
		for(int i=0; i<NLt; i++)
		{
			m_tPGsampler[i] = new PolyaGamma();
		}
		userNeighbour = myuserNeighbour;
		ind_fea = myindfea;
		Graphs g = new Graphs();
		userInNeighbours = g.generateuserInNeighbour(userNeighbour,U);
		tweetInNeighbours = g.generatetweetInNeighbour(users);
		retweet = new HashSet<String>();
		retweet = myretweet;
		cleanTempPrmts(users);
		computeTempPrmts(users, Z,C);
		computeSum(users, T);
		ptopic = new double[T];
		return true;
	}
	
	/**
	 * initialize the model from files
	 */
	protected boolean init_frIter(ArrayList<User> users, int vocabSize, int tsNo,HashMap<Integer, List<Integer>> myuserNeighbour, double[][] myindfea, String iterfile, HashSet<String> myretweet) {
		U = users.size();
		V = vocabSize;
		M = tsNo;
		// assign topics and communities randomly
		Z = new int[U][];
		C = new int[U][];
		for (int i = 0; i < U; i++) {
			Z[i] = new int[users.get(i).getDocWords().length];
			C[i] = new int[Z[i].length];
		}
		
		FileUtil.readMatrix(iterfile+"Z", Z);
		FileUtil.readMatrix(iterfile+"C", C);
		m_uPGsampler = new PolyaGamma[NLu];
		for(int i=0; i<NLu; i++)
		{
			m_uPGsampler[i] = new PolyaGamma();
		}
		m_tPGsampler = new PolyaGamma[NLt];
		for(int i=0; i<NLt; i++)
		{
			m_tPGsampler[i] = new PolyaGamma();
		}
		userNeighbour = myuserNeighbour;
		ind_fea = myindfea;
		Graphs g = new Graphs();
		userInNeighbours = g.generateuserInNeighbour(userNeighbour,U);
		tweetInNeighbours = g.generatetweetInNeighbour(users);
		retweet = new HashSet<String>();
		retweet = myretweet;
		cleanTempPrmts(users);
		//computeTempPrmts(users, Z,C);
		FileUtil.readMatrix(iterfile+"nuc", NUC);
		FileUtil.readMatrix(iterfile+"nct", NCT);
		FileUtil.readMatrix(iterfile+"ntw", NTW);
		FileUtil.readMatrix(iterfile+"popu", popularity);
		computeSum(users, T);
		
		FileUtil.readMatrix(iterfile+"theta", theta);
		FileUtil.readMatrix(iterfile+"vPhi", vPhi);
		FileUtil.readMatrix(iterfile+"nu", nu);
		FileUtil.readMatrix(iterfile+"pai", pai);
		FileUtil.readMatrix(iterfile+"eta", eta);
		FileUtil.readMatrix(iterfile+"lambda", lambda);
		FileUtil.readMatrix(iterfile+"delta", delta);
		
		return true;
	}

	//clean value for variational parameters
	public void cleanTempPrmts(ArrayList<User> users) {
		// initial parameters 
		NTW = new int[T][];
		vPhi = new float[T][];
		popularity = new int[T][];
		eta = new double[T][][];
		for (int t = 0; t < T; t++) {
			NTW[t] = new int[V];
			vPhi[t] = new float[V];
			for (int i = 0; i < V; i++) {
				NTW[t][i] = 0;
				vPhi[t][i] = 0.0f;
			}
			popularity[t] = new int[M];
			for (int i = 0; i < M; i++) {
				popularity[t][i] = 0;
			}
			eta[t] = new double[A][];
			for(int i=0; i<A;i++)
			{
				eta[t][i] = new double[A];
				for(int j=0; j<A; j++)
					eta[t][i][j] = 1.0f;
			}			
		}

		NCT = new int[A][];
		theta = new float[A][];
		for (int i = 0; i < A; i++) {
			NCT[i] = new int[T];
			theta[i] = new float[T];			
			for (int t = 0; t < T; t++) {
				NCT[i][t] = 0;
				theta[i][t] = 0.0f;
			}
		}
		
		NUC = new int[U][];
		pai = new float[U][];
		for (int i = 0; i < U; i++) {
			NUC[i] = new int[A];
			pai[i] = new float[A];
			for (int t = 0; t < A; t++) {
				NUC[i][t] = 0;
				pai[i][t] = 0.0f;
			}
		}
		nu = new double[N2f];
		for(int i=0; i<N2f; i++)
		{
			nu[i] = 1.0d;
		}
		lambda = new float[NLu];
		for(int i=0; i<NLu;i++)
		{
			lambda[i] = 1.0f;
		}
		delta = new float[NLt];
		for(int i=0; i<NLt;i++)
		{
			delta[i] = 1.0f;
		}
	}
	
	//clean value for variational parameters
	public void cleanCounter(ArrayList<User> users) {
		// initial parameters NW[] NWT[][] NIT[][] NY[] NUT[][]
		NTW = new int[T][];
		popularity = new int[T][];
		for (int t = 0; t < T; t++) {
			NTW[t] = new int[V];
			for (int i = 0; i < V; i++) {
				NTW[t][i] = 0;
			}
			popularity[t] = new int[M];
			for (int i = 0; i < M; i++) {
				popularity[t][i] = 0;
			}
		}

		NCT = new int[A][];
		for (int i = 0; i < A; i++) {
			NCT[i] = new int[T];
			for (int t = 0; t < T; t++) {
				NCT[i][t] = 0;
			}
		}
		
		NUC = new int[U][];
		for (int i = 0; i < U; i++) {
			NUC[i] = new int[A];
			for (int t = 0; t < A; t++) {
				NUC[i][t] = 0;
			}
		}
	}
	
	public void printnoretweet(ArrayList<User> users, String outputDir)
	{
		ArrayList<String> retweet = new ArrayList<String>();
		for (int u = 0; u < U; u++) {
			for (int n = 0; n < users.get(u).getDocWords().length; n++) {
				if(!getNoretweets(u, n, users))
					retweet.add(u+"\t"+n);
			}			
		}
		FileUtil.writeLines(outputDir+"retweet.txt", retweet);
	}
	
	public void inference(int iteration, ArrayList<User> users,
			int saveStep, int saveTimes, String outputDir, int initer) {
		if (iteration < saveStep * saveTimes) {
			System.err.println("iteration should be at least: " + saveStep * saveTimes);
			System.exit(0);
		}
		niters = iteration;
		double arho = A * rho, talpha = T * alpha, vbeta = V * beta;
		for (int i = 0; i < iteration; i++) {
			System.out.println("iteration " + i);
			long begintime = System.nanoTime();			
			drawEtaCOLD(users);
//			long ts1= System.nanoTime();
//			System.out.println("draw eta: "+(ts1-begintime));
			drawLambda(users);
//			long ts2= System.nanoTime();
//			System.out.println("draw lambda: "+(ts2-ts1));
			drawDelta(users);
//			long ts3= System.nanoTime();
//			System.out.println("draw delta: "+(ts3-ts2));
			
//			//multithreading
//			ThreadEta te = new ThreadEta(this, users);
//			ThreadLambda tl = new ThreadLambda(this, users);
//			ThreadDelta td = new ThreadDelta(this, users);
//			
//			te.start();
//			tl.start();
//			td.start();
//			
//			try {
//				te.join();
//				tl.join();
//				td.join();
//			} catch (InterruptedException e1) {
//				// TODO Auto-generated catch block
//				e1.printStackTrace();
//			}
			
			int uRefIx = 0, tRefIx =0;
			for (int u = 0; u < U; u++) {
//				long beg2 = System.nanoTime();
				double[] gVal = ComputeNeighborLhoodG(users, u, uRefIx);
//				long end2= System.nanoTime();
//				System.out.println("draw gVal: "+(end2-beg2));
				for (int n = 0; n < users.get(u).getDocWords().length; n++) {					
					//sample p(c_un=c)	
//					long beg3 = System.nanoTime();
					SampleCommunity(users.get(u)
							.getDocTimeStamp()[n], u, n, users, tRefIx, arho, talpha, gVal);
//					long end3= System.nanoTime();
//					System.out.println("draw community: "+(end3-beg3));
					//sample p(z_un=z)
//					long beg4 = System.nanoTime();
					SampleTopic(users.get(u).getDocWords()[n], users.get(u)
							.getDocTimeStamp()[n], u, n, users, tRefIx,talpha, vbeta);
//					long end4= System.nanoTime();
//					System.out.println("draw topic: "+(end4-beg4));
					if(users.get(u).getRetweets()!=null)
					{
						if(users.get(u).getRetweets().containsKey(n))
						{
							tRefIx++;
						}
					}
				}
				if(userNeighbour.get(u)!=null)
				{
					uRefIx += userNeighbour.get(u).size();
				}
			}	
			
			if (i >= (iteration - 1 - (saveStep * (saveTimes - 1))))
			{
				if ((iteration - i - 1) % saveStep == 0) {
					System.out.println("Saveing the mod el at " + (i + 1) + "-th iteration");
					uRefIx = 0; 
					tRefIx =0;
					for (int u = 0; u < U; u++) {
						//approximation, here we do not consider gVal for each single tweet, NUC[u][c]-- was need before the gVal calculation
						double[] gVal = ComputeNeighborLhoodG(users, u, uRefIx);

						for (int n = 0; n < users.get(u).getDocWords().length; n++) {
							//sample p(c_un=c)
							SampleCommunity_Final(users.get(u)
									.getDocTimeStamp()[n], u, n, users, tRefIx,arho, talpha,gVal);
							//sample p(z_un=z)
							SampleTopic_Final(users.get(u).getDocWords()[n], users.get(u)
									.getDocTimeStamp()[n], u, n, users, tRefIx,talpha, vbeta);
						}
						if(userNeighbour.get(u)!=null)
						{
							uRefIx += userNeighbour.get(u).size();
						}
						if(users.get(u).getRetweets()!=null)
							tRefIx += users.get(u).getRetweets().size();
					}
/*					List<int[]> instances = new ArrayList<int[]>();
					instances = getInstance("data/toy2/TMinput/");
					LearnNu(instances, users);*/
					try {
						computeModelParameter();
						saveModel(outputDir,initer+i + 1);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			long endtime = System.nanoTime();
			System.out.println(endtime-begintime);			
		}
		
		//PrintRetweetPro(users);
	}

	public void computeModelParameter() {
		System.out.println("computing model parameters...");

		for (int t = 0; t < T; t++) {
			for (int w = 0; w < V; w++)
				vPhi[t][w] = (float) ((NTW[t][w] + beta) / (SNTW[t] + V * beta));
		}

		for (int a = 0; a < A; a++) {
			for (int t = 0; t < T; t++) {
				theta[a][t] = (float) ((NCT[a][t] + alpha) / (SNCT[a] + T * alpha));
			}
		}
		
		for(int i=0; i<U;i++) {
			for(int a=0; a<A;a++) {
				pai[i][a] = (float) ((NUC[i][a] + rho) / (SNUC[i] + A*rho));
			}
		}		
		
		System.out.println("model parameters are computed");
	}

	private void computeSum(ArrayList<User> users, int T) {
		SNUC = new float[users.size()];
		for (int i = 0; i < users.size(); i++) {
			SNUC[i] = MatrixUtil.sumRow(NUC, i);
		}
		SNTW = new float[T];
		for (int t = 0; t < T; t++) {
			SNTW[t] = MatrixUtil.sumRow(NTW, t);
		}
		SNCT = new float[A];
		for(int a=0; a<A; a++)
		{
			SNCT[a]=MatrixUtil.sumRow(NCT, a);
		}
		SNTC = new float[T];
		for(int t=0; t<T; t++)
		{
			SNTC[t] = MatrixUtil.sumColumn(NCT, t, A);
		}
		SPTM = new float[T];
		for(int t=0; t<T; t++)
		{
			SPTM[t] = MatrixUtil.sumRow(popularity,t);
		}
	}

	private void computeTempPrmts(ArrayList<User> users, int[][] newZ, int[][] newC) {
		for (int i = 0; i < U; i++) {
			for (int j = 0; j < users.get(i).getDocWords().length; j++) {
				NUC[i][newC[i][j]]++;
				NCT[newC[i][j]][newZ[i][j]]++;
				for (int k = 0; k < users.get(i).getDocWords()[j].length; k++)
					NTW[newZ[i][j]][users.get(i).getDocWords()[j][k]]++;
				popularity[Z[i][j]][users.get(i).getDocTimeStamp()[j]]++;
			}
		}
	}

	public void getResFromLastIteration(ArrayList<User> users) {
		System.out.println("getting results from last interation...");
		//cleanTempPrmts(users);
		cleanCounter(users);
		computeTempPrmts(users, Z,C);
	}

//	private double[] ComputeNeighborLhoodG(ArrayList<User> users, int userInd, int RefIx)
//	{
//	//	long begintime = System.nanoTime();
//		double[] dLhoodVal = new double[A];
//		double sumval = 0D;
//		for(int a=0; a<A; a++)
//		{
//			dLhoodVal[a] = 0;
//			int community = a;
//			List<Integer> neighbours = new ArrayList<Integer>();
//			NUC[userInd][community] ++;
//			
//			if(userNeighbour.get(userInd)!=null)
//			{
//				neighbours = userNeighbour.get(userInd);	//neighbours of user userInd, the one who she follows
//				for(int i=0; i< neighbours.size(); i++)
//				{
//					int uNeighbor = neighbours.get(i);
//					//pai_u * pai_v
//					double weight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
//					dLhoodVal[a] += (weight - lambda[RefIx+i]*weight*weight);
//				}
//			}		
//			
//			if(userInNeighbours.get(userInd)!=null)
//			{
//				List<int[]> inneighbours = new ArrayList<int[]>();
//				//for all the neighbours of user userInd, the one who follow her
//				inneighbours = userInNeighbours.get(userInd);
//				for(int i=0; i< inneighbours.size(); i++)
//				{
//					int uNeighbor = inneighbours.get(i)[0];
//					//pai_u * pai_v
//					double weight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
//					dLhoodVal[a] += (weight - lambda[inneighbours.get(i)[1]]*weight*weight);
//				}	
//			}
//			dLhoodVal[a] = Math.exp(3*dLhoodVal[a]);
//			NUC[userInd][community] --;
//			sumval += dLhoodVal[a];
//		}
//		
//		for(int a=0; a<A; a++)
//		{
//			dLhoodVal[a] = dLhoodVal[a]/sumval;
//		}
////		long endtime = System.nanoTime();
////		System.out.println("compute G: "+(endtime-begintime));  
//		return dLhoodVal;
//	}

	private double[] ComputeNeighborLhoodG(ArrayList<User> users, int userInd, int RefIx)
	{
//		long begintime = System.nanoTime();
		double[] dLhoodVal = new double[A];
		double sumval = 0D;
		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = 0D;
		}

			List<Integer> neighbours = new ArrayList<Integer>();
			
			if(userNeighbour.get(userInd)!=null)
			{
				neighbours = userNeighbour.get(userInd);	//neighbours of user userInd, the one who she follows
				for(int i=0; i< neighbours.size(); i++)
				{
					int uNeighbor = neighbours.get(i);
					//pai_u * pai_v
					//double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
					double unitlength = (double)(users.get(userInd).getDocWords().length* users.get(uNeighbor).getDocWords().length);
					double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor])/unitlength;  // discriminant function value between pair <d, nNeighbor>
					for(int a=0; a<A; a++)
					{
						double weight = tweight+ (NUC[uNeighbor][a]/unitlength);
						//weight/= (double)(users.get(userInd).getDocWords().length* users.get(uNeighbor).getDocWords().length);
						dLhoodVal[a] += (weight - lambda[RefIx+i]*weight*weight);
					}
				}
			}		
			
			if(userInNeighbours.get(userInd)!=null)
			{
				List<int[]> inneighbours = new ArrayList<int[]>();
				//for all the neighbours of user userInd, the one who follow her
				inneighbours = userInNeighbours.get(userInd);
				for(int i=0; i< inneighbours.size(); i++)
				{
					int uNeighbor = inneighbours.get(i)[0];
					//pai_u * pai_v
					//double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor], users.get(userInd).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
					double unitlength = (double)(users.get(userInd).getDocWords().length* users.get(uNeighbor).getDocWords().length);
					double tweight = uDiscFun(NUC[userInd], NUC[uNeighbor])/ unitlength;  // discriminant function value between pair <d, nNeighbor>
					for(int a=0; a<A; a++)
					{
						double weight = tweight + (NUC[uNeighbor][a]/unitlength);
						//weight/= (double)(users.get(userInd).getDocWords().length* users.get(uNeighbor).getDocWords().length);
						dLhoodVal[a] += (weight - lambda[inneighbours.get(i)[1]]*weight*weight);
					}
				}	
			}
			
		//}
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = Math.exp(3*dLhoodVal[a]);	
			sumval += dLhoodVal[a];
		}
			
		//System.out.print("G:");
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = dLhoodVal[a]/sumval;
			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
			{
				for(int i=0; i<A; i++)
				{
					dLhoodVal[i] = 1;
				}
				return dLhoodVal;
			}
//			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
//				dLhoodVal[a] = 1;
		}

		//System.out.println("");
//		long endtime = System.nanoTime();
//		System.out.println("compute G: "+(endtime-begintime));  
		return dLhoodVal;
	}
	
	private double[] ComputeNeighborLhoodE_C(ArrayList<User> users, int userInd, int tweetInd, int RefIx, int topic)
	{
		//long begintime = System.nanoTime();
		double[] dLhoodVal = new double[A];
		HashMap<Integer, int[]> retweetinfo = users.get(userInd).getRetweets();
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{userInd, tweetInd};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = 0D;
		}
		double sumval = 0D;
			

		//for the tweet where this tweet retweet, if any
		if(retweetinfo!=null)
		{
			if(retweetinfo.containsKey(tweetInd))
			{
				int v = retweetinfo.get(tweetInd)[0];	//the user where the tweet is from
				int z = Z[userInd][tweetInd];	//the topic label of this retweet
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(userInd,v);
//					f=IndFea(v);
				//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
				float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z);
				double unitlength = (double) (users.get(userInd).getDocWords().length*users.get(v).getDocWords().length);
				for(int a=0; a<A; a++)
				{
					double weight = tweight;
					for (int cp=0; cp<A; cp++)
					{
						weight += eta[z][a][cp]*NUC[v][cp];
					}
					weight = (weight * para)/unitlength + popu*parapopu;
					dLhoodVal[a] += (weight - delta[RefIx]*weight*weight);
				}
			}
		}
		
		//for all tweets which retweet current tweet
		if(intweetneighbours!=null)
		{
			for(int i=0; i<intweetneighbours.size();i++)
			{
				int v = intweetneighbours.get(i)[0];
				int vtIx = intweetneighbours.get(i)[1];
				int z = Z[v][vtIx];
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(v,userInd);
//					f=IndFea(userInd);
				float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				//dLhoodVal[a] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z);
				double unitlength = (double) (users.get(userInd).getDocWords().length*users.get(v).getDocWords().length);
				for(int a=0; a<A; a++)
				{
					double weight = tweight;
					for (int cp=0; cp<A; cp++)
					{
						weight += eta[z][a][cp]*NUC[v][cp];
					}
					weight = (weight * para)/unitlength + popu*parapopu;
					dLhoodVal[a] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
				}
			}
		}


		
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = Math.exp((dLhoodVal[a]));
			sumval+=dLhoodVal[a];
		}
		
//		System.out.print("EC:");
		for(int a=0; a<A; a++)
		{
			dLhoodVal[a] = dLhoodVal[a]/sumval;
			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
			{
				for(int i=0; i<A; i++)
				{
					dLhoodVal[i] = 1;
				}
				return dLhoodVal;
			}
//			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
//				dLhoodVal[a] = 1;
		//	System.out.print(dLhoodVal[a]+" ");
		}
	//	System.out.println("");
//		long endtime = System.nanoTime();
//		System.out.println("compute ec: "+(endtime-begintime));
		
		return dLhoodVal;
	}
	
//	private double[] ComputeNeighborLhoodE_C(ArrayList<User> users, int userInd, int tweetInd, int RefIx, int topic)
//	{
//		//long begintime = System.nanoTime();
//		double[] dLhoodVal = new double[A];
//		HashMap<Integer, int[]> retweetinfo = users.get(userInd).getRetweets();
//		//for all tweets which retweet current tweet
//		List<int[]> intweetneighbours = new ArrayList<int[]>();
//		int[] upair = new int[]{userInd, tweetInd};
//		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
//		
//		
//		double sumval = 0D;
//		for(int a=0; a<A; a++)
//		{
//			dLhoodVal[a] = 0D;
//			int community = a;
//			NCT[community][topic] ++;
//			NUC[userInd][community]++;
//			//for the tweet where this tweet retweet, if any
//			if(retweetinfo!=null)
//			{
//				if(retweetinfo.containsKey(tweetInd))
//				{
//					int v = retweetinfo.get(tweetInd)[0];	//the user where the tweet is from
//					int z = Z[userInd][tweetInd];	//the topic label of this retweet
////					double[] f= new double[2*Nf];
////					//f = concatIndFea(userInd,v);
////					f=IndFea(v);
//					//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
//					float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
//					double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
//					dLhoodVal[a] += (weight - delta[RefIx]*weight*weight);
//				}
//			}
//			
//			//for all tweets which retweet current tweet
//			if(intweetneighbours!=null)
//			{
//				for(int i=0; i<intweetneighbours.size();i++)
//				{
//					int v = intweetneighbours.get(i)[0];
//					int vtIx = intweetneighbours.get(i)[1];
//					int z = Z[v][vtIx];
////					double[] f= new double[2*Nf];
////					//f = concatIndFea(v,userInd);
////					f=IndFea(userInd);
//					float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
//					double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
//					dLhoodVal[a] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
//				}
//			}
//			NCT[community][topic] --;
//			NUC[userInd][community]--;	
//		}
//		
//		for(int a=0; a<A; a++)
//		{
//			dLhoodVal[a] = Math.exp((dLhoodVal[a]));
//			sumval+=dLhoodVal[a];
//		}
//		
//		for(int a=0; a<A; a++)
//		{
//			dLhoodVal[a] = dLhoodVal[a]/sumval;
//			if(dLhoodVal[a]==0 || Double.isNaN(dLhoodVal[a]) ||Double.isInfinite(dLhoodVal[a]))
//				dLhoodVal[a] = 1;
//		}
//		
////		long endtime = System.nanoTime();
////		System.out.println("compute ec: "+(endtime-begintime));
//		
//		return dLhoodVal;
//	}
	

//	private double[] ComputeNeighborLhoodE_T(ArrayList<User> users, int userInd, int tweetInd, int RefIx, int community, int timestamp)
//	{
//		//long begintime = System.nanoTime();
//		double[] dLhoodVal = new double[T];
//		HashMap<Integer, int[]> retweetinfo = users.get(userInd).getRetweets();
//		//for all tweets which retweet current tweet
//		List<int[]> intweetneighbours = new ArrayList<int[]>();
//		int[] upair = new int[]{userInd, tweetInd};
//		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
//		
//		double sumval = 0D;			
//		for(int t=0; t<T; t++)
//		{
//			int topic =t;
//			dLhoodVal[t]=0D;
//			NCT[community][topic] ++;
//			NUC[userInd][community]++;
//	
//			popularity[topic][timestamp]++;
//			SPTM[topic]++;			
//				
//			//for the tweet where this tweet retweet, if any
//			if(retweetinfo!=null)
//			{
//				if(retweetinfo.containsKey(tweetInd))
//				{
//					int v = retweetinfo.get(tweetInd)[0];	//the user where the tweet is from
//					int z = Z[userInd][tweetInd];	//the topic label of this retweet
////					double[] f= new double[2*Nf];
////					//f = concatIndFea(userInd,v);
////					f=IndFea(v);
//					//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
//					float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
//					double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
//					dLhoodVal[t] += (weight - delta[RefIx]*weight*weight);
//				}
//			}
//			
//			//for all tweets which retweet current tweet
//			if(intweetneighbours!=null)
//			{
//				for(int i=0; i<intweetneighbours.size();i++)
//				{
//					int v = intweetneighbours.get(i)[0];
//					int vtIx = intweetneighbours.get(i)[1];
//					int z = Z[v][vtIx];
////					double[] f= new double[2*Nf];
////					//f = concatIndFea(v,userInd);
////					f=IndFea(userInd);
//					float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
//					double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
//					dLhoodVal[t] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
//				}
//			}
//			NCT[community][topic] --;
//			NUC[userInd][community]--;
//	
//			popularity[topic][timestamp]--;
//			SPTM[topic]--;
//		}
//		for(int t=0; t<T; t++)
//		{
//			dLhoodVal[t] = Math.exp((dLhoodVal[t]));
//			sumval+=dLhoodVal[t];
//		}
//		
//		for(int t=0; t<T; t++)
//		{
//			dLhoodVal[t] = dLhoodVal[t]/sumval;
//			if(dLhoodVal[t]==0 || Double.isNaN(dLhoodVal[t]) ||Double.isInfinite(dLhoodVal[t]))
//				dLhoodVal[t] = 1;
//		}
//
////		long endtime = System.nanoTime();
////		System.out.println("compute et: "+(endtime-begintime));
//		return dLhoodVal;
//	}
	
	/*private double[] ComputeNeighborLhoodE_T(ArrayList<User> users, int userInd, int tweetInd, int RefIx, int community, int timestamp)
	{
		//long begintime = System.nanoTime();
		double[] dLhoodVal = new double[T];
		HashMap<Integer, int[]> retweetinfo = users.get(userInd).getRetweets();
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{userInd, tweetInd};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		
		for(int t=0; t<T; t++)
		{
			dLhoodVal[t] = 0D;
		}
		
		double sumval = 0D;			
			
		//for the tweet where this tweet retweet, if any
		if(retweetinfo!=null)
		{
			if(retweetinfo.containsKey(tweetInd))
			{
				int v = retweetinfo.get(tweetInd)[0];	//the user where the tweet is from
				int z = Z[userInd][tweetInd];	//the topic label of this retweet
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(userInd,v);
//					f=IndFea(v);
				//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
				
				//double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				//dLhoodVal[t] += (weight - delta[RefIx]*weight*weight);
				double tweight = tDiscFun(NUC[userInd],NUC[v],z);
				double unitlength = (double) (users.get(userInd).getDocWords().length*users.get(v).getDocWords().length);
				for(int t=0; t<T; t++)
				{
					double weight = (tweight* para)/unitlength;
					float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
					weight = popu*parapopu;
					dLhoodVal[a] += (weight - delta[RefIx]*weight*weight);
				}
			}
		}
		
		//for all tweets which retweet current tweet
		if(intweetneighbours!=null)
		{
			for(int i=0; i<intweetneighbours.size();i++)
			{
				int v = intweetneighbours.get(i)[0];
				int vtIx = intweetneighbours.get(i)[1];
				int z = Z[v][vtIx];
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(v,userInd);
//					f=IndFea(userInd);
				float popu =  (float)popularity[z][users.get(userInd).getDocTimeStamp()[tweetInd]]/SPTM[z];
				double weight = tDiscFun(NUC[userInd],NUC[v],users.get(userInd).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
				dLhoodVal[t] += (weight - delta[intweetneighbours.get(i)[2]]*weight*weight);
			}
		}

		for(int t=0; t<T; t++)
		{
			dLhoodVal[t] = Math.exp((dLhoodVal[t]));
			sumval+=dLhoodVal[t];
		}
		
		for(int t=0; t<T; t++)
		{
			dLhoodVal[t] = dLhoodVal[t]/sumval;
			if(dLhoodVal[t]==0 || Double.isNaN(dLhoodVal[t]) ||Double.isInfinite(dLhoodVal[t]))
				dLhoodVal[t] = 1;
		}

//		long endtime = System.nanoTime();
//		System.out.println("compute et: "+(endtime-begintime));
		return dLhoodVal;
	}*/
	
	//to record the tweets which are not retweets and haven't been retweeted by others
	private boolean getNoretweets(int u, int n, ArrayList<User> users)
	{
		boolean res = false;
		HashMap<Integer, int[]> retweetinfo = users.get(u).getRetweets();
		//for all tweets which retweet current tweet
		List<int[]> intweetneighbours = new ArrayList<int[]>();
		int[] upair = new int[]{u, n};
		intweetneighbours = tweetInNeighbours.get(ComUtil.generateKey(upair));
		
		if(retweetinfo==null&&intweetneighbours==null)
		{
			res = true;
		}
		else if(retweetinfo!=null && intweetneighbours==null)
		{
			if(!retweetinfo.containsKey(n))
				res = true;
		}
		return res;
	}
	
	private boolean SampleCommunity(int ts, int u, int n, ArrayList<User> users, int tRefIx,double arho, double talpha, double[] gVal) {
		
		int topic = Z[u][n];
		int community = C[u][n];

		NUC[u][community]--;
		SNUC[u]--;
		NCT[community][topic]--;
		SNCT[community]--;
		SNTC[topic]--;

		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[A];
		double NUCsumRowU = SNUC[u];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_C(users, u, n, tRefIx,topic);
		}
//		long st = System.nanoTime();
		for (int a = 0; a < A; a++) {			
			double p1 = (double) (NUC[u][a] + rho) / (NUCsumRowU + arho);			
			double p2 = (double) (NCT[a][topic] + alpha) / (SNCT[a] + talpha);
			if(noret)
				pt[a] = p1 * p2 * gVal[a];
			else
				pt[a] = p1 * p2 * gVal[a] * eVal[a];
		}
//		long et = System.nanoTime();
//		System.out.println("in sample c: "+(et-st));
		//System.out.print("Sample community ");
		// cummulate multinomial parameters
		int sample = ComUtil.sample(pt, A);
		assert (sample >= 0 && sample < A) : "sample value error:" + sample;

		C[u][n] = sample;
		community = sample;

		NUC[u][community]++;
		SNUC[u]++;
		NCT[community][topic]++;
		SNCT[community]++;
		SNTC[topic]++;
		
		return true;
	}
	
	public void SampleOneTopic(int community, int i, double NCTsumRowC, double talpha, ArrayList<Integer> tempUniqueWords, ArrayList<Integer> tempCounts, double vbeta)
	{

		int wcount = 0;
		double p1 = (double) (NCT[community][i] + alpha) / (NCTsumRowC + talpha);
		double p2 = 1.0D;
		for (int w = 0; w < tempUniqueWords.size(); w++) {
			int tempvalue = NTW[i][tempUniqueWords.get(w)];
			// double sumRow = MatrixUtil.sumRow(NTW, i);
			double NTWsumRowT = SNTW[i];
			// checkEqual(sumRow, MatrixUtil.sumRow(NTW, i), "NTW");
			for (int numC = 0; numC < tempCounts.get(w); numC++) {
				p2 = p2 * ((double) (tempvalue + beta + numC) / ((double) NTWsumRowT
								+ vbeta + wcount));
				wcount++;
			}
		}
//		if(noret)
			ptopic[i] = p1 * p2;
//		else
//		{
//			pt[i] = p1 * p2  * eVal[i];
//		}
	}
	
	private boolean SampleTopic(int[] words, int ts, int u, int n, ArrayList<User> users, int tRefIx, double talpha, double vbeta) {
		
		int topic = Z[u][n];
		int community = C[u][n];
		int timestamp = users.get(u).getDocTimeStamp()[n];
		// get words and their count in [u,n]
		ArrayList<Integer> tempUniqueWords = new ArrayList<Integer>();
		ArrayList<Integer> tempCounts = new ArrayList<Integer>();
		uniqe(words, tempUniqueWords, tempCounts);

		NCT[community][topic]--;
		SNCT[community]--;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] -= tempCounts.get(w1);
			SNTW[topic] -= tempCounts.get(w1);
		}
		SNTC[topic]--;
		popularity[topic][timestamp] --;
		SPTM[topic]--;
		SNTC[topic]--;
		
		// get p(Z_{u,n} = z|Z_c, W, Y, I)
//		double[] pt = new double[T];
		double NCTsumRowC = SNCT[community];
//		boolean noret = false;
//		double[] eVal = null;
//		if(noretweet.contains(u+"\t"+n))
//		{
//			noret = true;
//		}
//		else
//		{			
//			eVal = ComputeNeighborLhoodE_T(users, u, n, tRefIx,community,timestamp);
//		}
		
//		long st = System.nanoTime();
		for (int i = 0; i < T; i++) {
			int wcount = 0;
			double p1 = (double) (NCT[community][i] + alpha) / (NCTsumRowC + talpha);
			double p2 = 1.0D;
			for (int w = 0; w < tempUniqueWords.size(); w++) {
				int tempvalue = NTW[i][tempUniqueWords.get(w)];
				// double sumRow = MatrixUtil.sumRow(NTW, i);
				double NTWsumRowT = SNTW[i];
				// checkEqual(sumRow, MatrixUtil.sumRow(NTW, i), "NTW");
				for (int numC = 0; numC < tempCounts.get(w); numC++) {
					p2 = p2 * ((double) (tempvalue + beta + numC) / ((double) NTWsumRowT
									+ vbeta + wcount));
					wcount++;
				}
			}
//			if(noret)
				ptopic[i] = p1 * p2;
//			else
//			{
//				pt[i] = p1 * p2  * eVal[i];
//			}
		}
		
//		//multithreading
//		ThreadTopic[] tt = new ThreadTopic[T];
//		for(int i=0; i<T; i++)
//		{
//			tt[i] = new ThreadTopic(this, community, i, NCTsumRowC, talpha, tempUniqueWords, tempCounts, vbeta); 
//			tt[i].start();
//			try {
//				tt[i].join();
//			} catch (InterruptedException e1) {
//				// TODO Auto-generated catch block
//				e1.printStackTrace();
//			}
//		}
		
		
//		ThreadTopic tt0 = new ThreadTopic(this, community, 0, NCTsumRowC, talpha, tempUniqueWords, tempCounts, vbeta);
//		ThreadTopic tt1 = new ThreadTopic(this, community, 1, NCTsumRowC, talpha, tempUniqueWords, tempCounts, vbeta);
//		ThreadTopic tt2 = new ThreadTopic(this, community, 2, NCTsumRowC, talpha, tempUniqueWords, tempCounts, vbeta);
//		tt0.start();
//		tt1.start();
//		tt2.start();
//		
//		try {
//			tt0.join();
//			tt1.join();
//			tt2.join();
//		} catch (InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

//		long et = System.nanoTime();
//		System.out.println("in sample t: "+(et-st));

		// cummulate multinomial parameters
		int sample = ComUtil.sample(ptopic, T);
		assert (sample >= 0 && sample < T) : "sample value error:" + sample;

		Z[u][n] = sample;
		topic = sample;

		// update NTW[T][W](y=1) NTI[T][M] NUT[U][T] in {u,n}
		NCT[community][topic]++;
		SNCT[community]++;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] += tempCounts.get(w1);
			SNTW[topic] += tempCounts.get(w1);
		}
		SNTC[topic]++;
		popularity[topic][timestamp] ++;
		SPTM[topic]++;
		SNTC[topic]++;
		tempUniqueWords.clear();
		tempCounts.clear();
		
		return true;
	}
	
	private boolean SampleCommunity_Final(int ts, int u, int n, ArrayList<User> users, int tRefIx, double arho, double talpha, double[] gVal) {
		int topic = Z[u][n];
		int community = C[u][n];
		
		NUC[u][community]--;
		SNUC[u]--;
		NCT[community][topic]--;
		SNCT[community]--;
		SNTC[topic]--;

		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[A];
		double NUCsumRowU = SNUC[u];
		boolean noret = false;
		double[] eVal = null;
		if(!retweet.contains(u+"\t"+n))
		{
			noret = true;
		}
		else
		{			
			eVal = ComputeNeighborLhoodE_C(users, u, n, tRefIx,topic);
		}
		
		for (int a = 0; a < A; a++) {
			double p1 = (double) (NUC[u][a] + rho) / (NUCsumRowU + arho);
			double p2 = (double) (NCT[a][topic] + alpha) / (SNCT[a]+ talpha);
			if(noret)
				pt[a] = p1 * p2 * gVal[a];
			else
				pt[a] = p1 * p2 * gVal[a] * eVal[a];
		}
		// cummulate multinomial parameters
		//int sample = ComUtil.sample(pt, A);
		int sample = ComUtil.findMax(pt, A);
		assert (sample >= 0 && sample < A) : "sample value error:" + sample;

		C[u][n] = sample;
		community = sample;
		NUC[u][community]++;
		SNUC[u]++;
		NCT[community][topic]++;
		SNCT[community]++;
		SNTC[topic]++;
		return true;
	}
	
	private boolean SampleTopic_Final(int[] words, int ts, int u, int n, ArrayList<User> users, int tRefIx,double talpha, double vbeta) {
		int topic = Z[u][n];
		int community = C[u][n];
		int timestamp = users.get(u).getDocTimeStamp()[n];
		// get words and their count in [u,n]
		ArrayList<Integer> tempUniqueWords = new ArrayList<Integer>();
		ArrayList<Integer> tempCounts = new ArrayList<Integer>();
		uniqe(words, tempUniqueWords, tempCounts);

		NCT[community][topic]--;
		SNCT[community]--;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] -= tempCounts.get(w1);
			SNTW[topic] -= tempCounts.get(w1);
		}
		SNTC[topic]--;
		popularity[topic][timestamp] --;
		SPTM[topic]--;
		SNTC[topic]--;
		
		// get p(Z_{u,n} = z|Z_c, W, Y, I)
		double[] pt = new double[T];
		double NCTsumRowC = SNCT[community];
//		boolean noret = false;
//		double[] eVal = null;
//		if(noretweet.contains(u+"\t"+n))
//		{
//			noret = true;
//		}
//		else
//		{			
//			eVal = ComputeNeighborLhoodE_T(users, u, n, tRefIx,community,timestamp);
//		}
		for (int i = 0; i < T; i++) {
			int wcount = 0;
			double p1 = (double) (NCT[community][i] + alpha) / (NCTsumRowC + talpha);
			double p2 = 1.0D;
			for (int w = 0; w < tempUniqueWords.size(); w++) {
				int tempvalue = NTW[i][tempUniqueWords.get(w)];
				double NTWsumRowT = SNTW[i];
				for (int numC = 0; numC < tempCounts.get(w); numC++) {
					p2 = p2 * ((double) (tempvalue + beta + numC) / ((double) NTWsumRowT
									+ vbeta + wcount));
					wcount++;
				}
			}
//			if(noret)
				pt[i] = p1 * p2;
//			else
//				pt[i] = p1 * p2  * eVal[i];
		}

		// cummulate multinomial parameters
		//int sample = ComUtil.sample(pt, T);
		int sample = ComUtil.findMax(pt,T);
		assert (sample >= 0 && sample < T) : "sample value error:" + sample;

		Z[u][n] = sample;
		topic = sample;

		// update NTW[T][W](y=1) NTI[T][M] NUT[U][T] in {u,n}
		NCT[community][topic]++;
		SNCT[community]++;
		for (int w1 = 0; w1 < tempUniqueWords.size(); w1++) {
			NTW[topic][tempUniqueWords.get(w1)] += tempCounts.get(w1);
			SNTW[topic] += tempCounts.get(w1);
		}
		SNTC[topic]++;
		popularity[topic][timestamp] ++;
		SPTM[topic]++;
		SNTC[topic]++;
		tempUniqueWords.clear();
		tempCounts.clear();
		return true;
	}
	
	//lambda for user-user following friendship, assign lambda value for each user-user pair
	private void drawLambda(ArrayList<User> users)
	{
		double uDiscFuncVal;
		int uPIx = 0;
		for ( int u = 0; u < U; u++) {
			if(userNeighbour.get(u)!=null)
			{
				for ( int v=0; v<userNeighbour.get(u).size(); v++ ) {
					int uNeighbor = userNeighbour.get(u).get(v);
					//pai_u * pai_v
					uDiscFuncVal = uDiscFun(NUC[u], NUC[uNeighbor], users.get(u).getDocWords().length, users.get(uNeighbor).getDocWords().length);  // discriminant function value between pair <d, nNeighbor>
					lambda[uPIx] = (float) m_uPGsampler[uPIx].nextPG(1, uDiscFuncVal);
					uPIx ++;
				}
			}
		}
	}
	
	//delta for retweet, assign delta for each retweet
	private void drawDelta(ArrayList<User> users)
	{
		double tDiscFuncVal;
		int tPIx = 0;
		for(int u=0; u<U; u++)
		{
			HashMap<Integer, int[]> retweetinfo = users.get(u).getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keySet())
				{
					int v = retweetinfo.get(utIx)[0];	//the user where the tweet is from
					int z = Z[u][utIx];	//the topic label of this retweet
					
//					double[] f= new double[2*Nf];
//					//f = concatIndFea(u,v);
//					f=IndFea(v);
					float popu =  (float)popularity[z][users.get(u).getDocTimeStamp()[utIx]]/SPTM[z];
					tDiscFuncVal = tDiscFun(NUC[u],NUC[v],users.get(u).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu);
					if(!Double.isInfinite(tDiscFuncVal)&&!Double.isNaN(tDiscFuncVal))
					{
						delta[tPIx] = (float) m_tPGsampler[tPIx].nextPG(1, tDiscFuncVal);
					}
					tPIx ++;
				}
			}
		}
	}
	
	//delta for retweet, assign delta for each retweet
	private void drawEtaCOLD(ArrayList<User> users)
	{

		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					eta[t][i][j] = 0;
			}			
		}
		for(int u=0; u<U; u++)
		{
			HashMap<Integer, int[]> retweetinfo = users.get(u).getRetweets();
			if(retweetinfo!=null)
			{
				for(int utIx: retweetinfo.keySet())
				{
					int v = retweetinfo.get(utIx)[0];	//the user where the tweet is from
					int vtIx = retweetinfo.get(utIx)[1];
					int z = Z[u][utIx];	//the topic label of this retweet
					eta[z][C[u][utIx]][C[v][vtIx]] ++;					
				}
			}
		}
		
		//normalize to make eta<1
		double maxeta = 0;
		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					
				{
					if(eta[t][i][j]>maxeta)
						maxeta = eta[t][i][j];
				}
			}			
		}
		
		for (int t = 0; t < T; t++) {
			for(int i=0; i<A;i++)
			{
				for(int j=0; j<A; j++)
					eta[t][i][j] = eta[t][i][j]/maxeta;
			}			
		}
	}
	
	//concat the feature for user u and user v
	private double[] concatIndFea(int u, int v)
	{
		double[] f= new double[2*Nf];
		for(int i=0; i<Nf; i++)
		{
			f[i] = ind_fea[u][i]*paranu;
			f[Nf+i] = ind_fea[v][i]*paranu;
		}
		return f;
	}
	
	//use only the statistics of users being retweeted
	private double[] IndFea(int u)
	{
		double[] f= new double[Nf];
		for(int i=0; i<Nf; i++)
		{
			f[i] = ind_fea[u][i]*paranu;
		}
		return f;
	}
	
	
	//calculate c_hat u * c_hat v
	private double uDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v)
	{
		double sum = 0D;
		for(int i=0; i< c_hat_u.length; i++)
		{
			sum += (double) c_hat_u[i]*c_hat_v[i];
		}
		sum /= (double)(N_u*N_v);
		return sum;
	}
	
	//calculate c_hat u * c_hat v (nuc(u,c)++), this is to speed up ComputeNeighborLhoodG
	private double uDiscFun(int[] c_hat_u, int[] c_hat_v)
	{
		double sum = 0D;
		for(int i=0; i< c_hat_u.length; i++)
		{
			sum += (double) c_hat_u[i]*c_hat_v[i];
		}
		return sum;
	}
	
	//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f)
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v, float N_k, int topick, int[] z_hat_k, float popularity/*, double[] nu, double[] f*/)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]/**z_hat_k[c]*eta[c][cp]*/*eta[topick][c][cp]*c_hat_v[cp]/**z_hat_k[cp]*/);
			}
		}
		res *= para;
		res /= (double)(N_u*N_v/**N_k*N_k*/);
		res += popularity*parapopu;
		return res;
	}
	
	//calculate sum_c sum_c' (c_hat{u,c}z_hat{c,k}eta{c,c',k}c_hat{v,c'}z_hat{c',k}+popularity{k,t}+nu*f), to speed up computeNeighbourhoodEC
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int topick)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]*eta[topick][c][cp]*c_hat_v[cp]);
			}
		}
		return res;
	}
	
	private double tDiscFun(int[] c_hat_u, int[] c_hat_v, int N_u, int N_v, float N_k, int topick, int[] z_hat_k, float popularity, double[] nu, double[] f)
	{
		double res = 0D;
		for(int c=0; c<A; c++)
		{
			for(int cp=0; cp<A; cp++)
			{
				res+=(double)(c_hat_u[c]/**z_hat_k[c]*eta[c][cp]*/*eta[topick][c][cp]*c_hat_v[cp]/**z_hat_k[cp]*/);
			}
		}
		res *= para;
		res /= (double)(N_u*N_v/**N_k*N_k*/);
		res += (popularity*parapopu+MatrixUtil.vectorTimes(nu, f));			//need to normalize to make the three factors in the same scale
		return res;
	}
	
	
	public static void uniqe(int[] words, 
			ArrayList<Integer> tempUniqueWords, ArrayList<Integer> tempCounts) {
		for (int i = 0; i < words.length; i++) {
			if (tempUniqueWords.contains(words[i])) {
				int index = tempUniqueWords.indexOf(words[i]);
				tempCounts.set(index, tempCounts.get(index) + 1);
			} else {
				tempUniqueWords.add(words[i]);
				tempCounts.add(1);
			}

		}
	}
	
	public double classify(int z, int u, int v, int t, double[] f, ArrayList<User> users)
	{
		double res = 0;
		float popu =  (float)popularity[z][t]/SPTM[z];
		double prob = tDiscFun(NUC[u],NUC[v],users.get(u).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu,nu,f);
		
		res = 1/(1+Math.exp(-prob));
		return res;
	}
	
	List<int[]> getInstance(String dir)
	{
		List<int[]> instances = new ArrayList<int[]>();
		ArrayList<String> Pinstance = new ArrayList<String>();
		ArrayList<String> Ninstance = new ArrayList<String>();
		FileUtil.readLines(dir+"PInstance.txt", Pinstance);
		FileUtil.readLines(dir+"NInstance.txt", Ninstance);
		for (String ins: Pinstance)
		{
			String[] terms = ins.split(",");
			int[] in = new int[5];
			for(int i=0; i<4; i++)
			{
				in[i] = Integer.parseInt(terms[i]);
			}
			in[4] = 1;
			instances.add(in);
		}
		for (String ins: Ninstance)
		{
			String[] terms = ins.split(",");
			int[] in = new int[5];
			for(int i=0; i<4; i++)
			{
				in[i] = Integer.parseInt(terms[i]);
			}
			in[4] = 0;
			instances.add(in);
		}
		
		return instances;
	}
	//logistic regression to learn \nu
	public void LearnNu(List<int[]> instances, ArrayList<User> users)
	{
		double rate = 0.00035;
		for (int n = 0; n < 2000; n++)
		{
			double lik = 0.0;
			for(int i=0; i<instances.size(); i++)
			{
				int[] ins = instances.get(i);
				double[] f= new double[2*Nf];
				//f=IndFea(ins[2]);
				f = concatIndFea(ins[0],ins[2]);
				double predicted = classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3], f,users);
				int label = ins[4];
				
				for (int j=0; j<nu.length; j++) {
                    nu[j] = nu[j] + rate * (label - predicted) * f[j];
                }
				 lik += label * Math.log(classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3], f,users)) + (1-label) * Math.log(1- classify(Z[ins[0]][ins[1]],ins[0],ins[2],ins[3],f, users));
			}
           // System.out.println("iteration: " + n + " " + Arrays.toString(nu) + " mle: " + lik);
		}
	}
	
	//predict if the retweet will happen
	public double ProbRetweet(int z, int u, int v, int t, ArrayList<User> users)
	{
		double res = 0;
		double[] f= new double[2*Nf];
		f = concatIndFea(u,v);
		//f=IndFea(v);
		float popu =  (float)popularity[z][t]/SPTM[z];
		double prob = tDiscFun(NUC[u],NUC[v],users.get(u).getDocWords().length, users.get(v).getDocWords().length,SNTC[z],z,MatrixUtil.getColumn(NCT, z),popu,nu,f);
		
		res = 1/(1+Math.exp(-prob));
		return res;
	}
	
	public void PrintRetweetPro(ArrayList<User> users)
	{
		for (int u = 0; u < U; u++) {	
			for (int v=0; v<U;v++) {
				if(u!=v)
				{
					for(int z=0; z<T;z++)
					{
						//for(int t=0; t<M;t++)
						//{
							System.out.println(/*"The probability of user "+(u+1)+" retweet user "+(v+1)+" on topic "+z+" is "+*/ProbRetweet(z, u, v, 1, users));
						//}
					}
					System.out.println();
				}
			}
		}
	}

	/**
	 * output Model paramters
	 * */
	public void outputModelRes() {
		// output Z
		System.out.println("Z[u][n]: ");
		MatrixUtil.printArray(Z);
	}

	public void outTaggedDoc(ArrayList<User> users,
			/*ArrayList<String> uniWordMap,*/ String outputDir) {
		ArrayList<String> datalines = new ArrayList<String>();
		for (int i = 0; i < users.size(); i++) {
			for (int j = 0; j < users.get(i).getDocWords().length; j++) {
				String tmpline = "Community "+C[i][j]+", Topic " + Z[i][j] + ": ";
				tmpline += users.get(i).getDocTimeStamp()[j] +" ";
				for (int k1 = 0; k1 < users.get(i).getDocWords()[j].length; k1++) {
						tmpline += /*uniWordMap.get(*/
								users.get(i).getDocWords()[j][k1]/*)*/ +" ";
				}
				datalines.add(tmpline);
				System.out.println(tmpline);
			}
			System.out.println("");
			FileUtil.writeLines(outputDir + users.get(i).getId(),
					datalines);
			datalines.clear();
		}
	}

	void saveModelRes(String string) throws Exception {
		BufferedWriter writer = null;
		writer = new BufferedWriter(new FileWriter(new File(string)));
		writer.write("Z[u][n]: \n");
		for (int i = 0; i < Z.length; i++) {
			for (int j = 0; j < Z[i].length; j++)
				writer.write(Z[i][j] + "\t");
			writer.write("\n");
		}
		writer.flush();
		writer.close();
	}

	public boolean saveModel(String output, int iter) throws Exception {
		output = output+"iter"+iter+"/";
		(new File(output)).mkdirs();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.theta")));
		ModelComFunc.writeData(theta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.pai")));
		ModelComFunc.writeData(pai, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.eta")));
		ModelComFunc.writeData(eta, writer);
		writer.close();

		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.vPhi")));
		ModelComFunc.writeData(vPhi, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nu")));
		ModelComFunc.writeData(nu, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.lambda")));
		ModelComFunc.writeData(lambda, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.delta")));
		ModelComFunc.writeData(delta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nuc")));
		ModelComFunc.writeData(NUC, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nct")));
		ModelComFunc.writeData(NCT, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.ntw")));
		ModelComFunc.writeData(NTW, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.popu")));
		ModelComFunc.writeData(popularity, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.Z")));
		ModelComFunc.writeData(Z, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.C")));
		ModelComFunc.writeData(C, writer);
		writer.close();

		return true;
	}
	
	public boolean saveModel(String output/*, ArrayList<String> uniWordMap*/) throws Exception {
		ArrayList<Integer> rankList = new ArrayList<Integer>();
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.theta")));
		ModelComFunc.writeData(theta, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.pai")));
		ModelComFunc.writeData(pai, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.eta")));
		ModelComFunc.writeData(eta, writer);
		writer.close();

		/*writer = new BufferedWriter(new FileWriter(new File(output
				+ "model-topic-words.txt")));
		for (int t = 0; t < vPhi.length; t++) {
			ComUtil.getTop(vPhi[t], rankList, 20);
			writer.write("Topic " + t + "\n");
			ModelComFunc.writeData(vPhi[t], uniWordMap, rankList, writer, "\t");
			rankList.clear();
		}
		writer.close();*/
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.vPhi")));
		ModelComFunc.writeData(vPhi, writer);
		writer.close();
		
		writer = new BufferedWriter(new FileWriter(new File(output
				+ "model.nu")));
		ModelComFunc.writeData(nu, writer);
		writer.close();

		return true;
	}


	public void output(ArrayList<User> users, ArrayList<String> uniWords,
			ArrayList<String> uniItems, String outputDir) {
		ArrayList<String> datalines = new ArrayList<String>();
		for (int i = 0; i < users.size(); i++) {
			for (int j = 0; j < users.get(i).getDocWords().length; j++) {
				String tmpline = "";
				for (int k1 = 0; k1 < users.get(i).getDocWords()[j].length; k1++) {
					tmpline += uniWords.get(users.get(i).getDocWords()[j][k1])
							+ " ";
				}
				tmpline += uniItems.get(users.get(i).getDocTimeStamp()[j]);
				datalines.add(tmpline);
			}
			FileUtil.writeLines(outputDir + users.get(i).getId() + ".txt",
					datalines);
			datalines.clear();
		}
	}
	
	class ThreadEta extends Thread{
		private Modelparatopic m;
		private ArrayList<User> us;
		ThreadEta(Modelparatopic m, ArrayList<User> us){
			this.m = m;
			this.us = us;
		}
		
		public void run()
		{
			m.drawEtaCOLD(us);
		}
	}
	
	class ThreadLambda extends Thread{
		private Modelparatopic m;
		private ArrayList<User> us;
		ThreadLambda(Modelparatopic m, ArrayList<User> us){
			this.m = m;
			this.us = us;
		}
		
		public void run()
		{
			m.drawLambda(us);
		}
	}
	
	class ThreadDelta extends Thread{
		private Modelparatopic m;
		private ArrayList<User> us;
		ThreadDelta(Modelparatopic m, ArrayList<User> us){
			this.m = m;
			this.us = us;
		}
		
		public void run()
		{
			m.drawDelta(us);
		}
	}
	
	class ThreadTopic extends Thread{
		private Modelparatopic m;
		int community;
		int i;
		double NCTsumRowC;
		double talpha;
		ArrayList<Integer> tempUniqueWords;
		ArrayList<Integer> tempCounts;
		double vbeta;
		
		ThreadTopic(Modelparatopic m, int community, int i, double NCTsumRowC, double talpha, ArrayList<Integer> tempUniqueWords, ArrayList<Integer> tempCounts, double vbeta){
			this.m = m;
			this.community = community;
			this.i = i;
			this.NCTsumRowC = NCTsumRowC;
			this.talpha = talpha;
			this.tempUniqueWords = tempUniqueWords;
			this.tempCounts = tempCounts;
			this.vbeta = vbeta;
		}
		
		public void run()
		{
			m.SampleOneTopic(community, i, NCTsumRowC, talpha, tempUniqueWords, tempCounts, vbeta);
		}
	}

}
