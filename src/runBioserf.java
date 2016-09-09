import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Collections;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ByteArrayInputStream;
import java.io.UnsupportedEncodingException;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.channels.FileChannel;
import org.xml.sax.SAXException;
import org.xml.sax.InputSource;
import org.biojava.bio.program.sax.BlastLikeSAXParser;
import org.biojava.bio.program.ssbind.BlastLikeSearchBuilder;
import org.biojava.bio.program.ssbind.SeqSimilarityAdapter;
import org.biojava.bio.search.*;
import org.biojava.bio.seq.db.*;
import org.biojava.bio.symbol.Alignment;
import org.ucl.util.Utils;
import org.ucl.util.ExternalProcess;

class runBioserf
{
  static String strJobExt = "bioserf";
  static String seq;

  public static void main( String args[] ) throws Exception
  {
    if( args.length != 17 )
    {
          System.out.println( "usage: java -cp bioserf/src/:bioserf/src/org/ucl/ parseBlast blast_file seq_file evalueThreshold(0.0005) /tmp/ PDB_PATH PDBAA_PATH percentIDThreshold(40) strUUID mod9.17 B0R5N0.pgen.presults /modeller9.17/modlib/ /modeller9.17/lib/x86_64-intel8/ /usr/bin/python qmodcheck_mainens qmodcheck modcheckpot.dat tmjury3d_mq_modeller" );
          System.exit(1);
    }

    String strBlastFile = args[0]; //the input blast file, blast vs pdbaa
    String strSeqFile = args[1];   //the blastquery sequence
    Double eValueThreshold = Double.parseDouble(args[2]); //blast eValueThreshold for selection
    String strTmpPath = args[3]; //path to write output files
    String strPDB = args[4];     //location of all your pdb files
    String strPDBAA = args[5];   //location of the pdbaa
    Double percentIDThreshold = Double.parseDouble(args[6]); //blast percent ID for selection
    String strUUID = args[7];   //an id to uniquely id the files
    String strModellerBin = args[8];
    String strPResults = args[9];
    String strModLib = args[10];
    String strModArch = args[11];
    String strPython = args[12];
    String strQModCheck_mainens = args[13];
    String strQModCheck = args[14];
    String strModCheckPotDat = args[15];
    String strTmJury = args[16];

    //read in blast data
    File fBlastOut = new File(strBlastFile);
    File fPResults = new File(strPResults);
    File fFinalModel = new File(strTmpPath + strUUID + ".B99990001.pdb");

    seq = readSequence(strSeqFile);
    String strBlastData = Utils.fileToStr(fBlastOut.getCanonicalPath());
    outputPsiblastModellerInput(strBlastData, eValueThreshold,
                                strTmpPath, strPDB, strPDBAA,
                                strUUID, percentIDThreshold);
    runModeller(strUUID, "psiblast", strTmpPath, strModellerBin);
    System.out.println("Built PSI-BLAST models!\n");
    System.out.println("Tidying pGenTHREADER models!\n");
    tidyGenModels(fPResults, strTmpPath, strUUID);

    File fEnsemble = catModels(strTmpPath, strUUID);

    //Run David's model collation script, read in the final model and push it to the frontend
    //TODO: Add chain ID to the output pdb file
    buildFinalModel(strTmpPath, strUUID, strModLib, strModArch, strPython,
                    fEnsemble, strSeqFile, strQModCheck_mainens, strQModCheck,
                    strModCheckPotDat, strTmJury);

    //go to the tmp directory
    String strFinalModel = Utils.fileToStr(fFinalModel.getCanonicalPath());
    strFinalModel = toCASPFormat(strFinalModel, strUUID);
    Utils.strToFile(strUUID+".final.pdb", strFinalModel);
  }

  //TODO: add a chain ID if blank, for resubmission tasks
   private static String toCASPFormat(String strFinalModel, String strSeqName) {

      String[] lines = strFinalModel.split("\n");

      // Get the target name

      String target_name = strSeqName;
      String author_code = "NewSerf";
      String remark_text = "";

      // Trim out what we don't need

      StringBuilder sb = new StringBuilder();

      sb.append("PFRMAT TS" + "\n");
      sb.append("TARGET " + target_name + "\n");
      sb.append("AUTHOR " + author_code + "\n");
      sb.append("REMARK " + remark_text + "\n");
      sb.append("PARENT N/A\n");
      sb.append("METHOD   1  PSIBLAST vs to find close homologos\n");
      sb.append("METHOD   2  pGenTHREADER to recognize distant homologs/folds\n");
      sb.append("METHOD   3  HHBlitz to recognize folds distant homologs\n");
      sb.append("METHOD   4  all-against-all TM score, select the maximally connected homologs\n");
      sb.append("METHOD   5  Uses " + "\n");
      sb.append("MODEL    1\n");

      for (int i = 0; i < lines.length; i++) {
          if (lines[i].startsWith("ATOM") || lines[i].startsWith("HETATM")) {
              // Necessary to edit Modeller's occupancy
              // field since it's not always 1.0
              char[] chararr = lines[i].toCharArray();

              // Offset of occupancy

              int occoff = 56;
              String template_occ = "1.00  0.00";

              for (int k = 0; k < template_occ.length(); k++) {
                  chararr[occoff + k] = template_occ.charAt(k);
              }

              String correctedline = new String(chararr);
              sb.append(correctedline + "\n");
          }
      }
      sb.append("TER" + "\n");
      sb.append("END" + "\n");

      return sb.toString();
  }

  protected static void RunQModCheck(String strQModCheck_mainens, String strQModCheck, String strModCheckPotDat, String strTmpPath, String strEnsemble, String strQmo) throws Exception
  {
      try
      {
          //../bin/qmodcheck_mainens B0R5N0.ensemble.pdb ../bin/qmodcheck ../data/modcheckpot.dat ~/Code/bioserf/src/
          String strCurCmd = strQModCheck_mainens + " " + strEnsemble +" "+strQModCheck+" "+strModCheckPotDat+" "+strTmpPath;
          ExternalProcess epQM = new ExternalProcess(strCurCmd.split(" "));
          System.out.println("    Running qmodcheck_mainens\n    " + strCurCmd);
          int nRes = epQM.call();
          Thread.sleep(1000);
          if (nRes != 0) {
              //job.UpdateStatus(StatusClass.Error, StatusCode.JobRunning, "psipass2 failed", strError);
              throw new Exception("Error calling qmodcheck_mainens");
          }
          String strRes = epQM.getOutput();
          Utils.strToFile(strQmo, strRes);
      }
      catch (Exception ex)
      {
          //job.UpdateStatus(StatusClass.Error, StatusCode.JobRunning, "Psipass call failed", strError);
          throw (new Exception(ex.getMessage()));
      }
  }
   protected static void RunTmJury(String strTmJury, String strTmp, String strPrefix, String strFasta, String strQmo, String strEnsemble, String strModelCount) throws Exception
  {
      // /var/www/cgi-bin/psipred/bin/tmjury3d_mq_modeller $tempdir prefix $tempdir/$jobid.fasta $tempdir/$jobid.mg3d.qmo $tempdir/$jobid.mg3d.ensemble.pdb 10 > /dev/null") != 0)
      try
      {
          String strCurCmd = strTmJury + " " + strTmp + " "+ strPrefix + " " + strFasta + " " + strQmo + " "  + strEnsemble + " " + strModelCount;
          ExternalProcess epTJ = new ExternalProcess(strCurCmd.split(" "));
          System.out.println("    Running tmjury\n    " + strCurCmd);
          int nRes = epTJ.call();
          Thread.sleep(1000);
          if (nRes != 0) {
              //job.UpdateStatus(StatusClass.Error, StatusCode.JobRunning, "psipass2 failed", strError);
              throw new Exception("Error calling tmjury");
          }
      }
      catch (Exception ex)
      {
          //job.UpdateStatus(StatusClass.Error, StatusCode.JobRunning, "Psipass call failed", strError);
          throw (new Exception(ex.getMessage()));
      }
  }
  protected static void buildFinalModel(String strTmpPath, String strUUID, String strModLib, String strModArch, String strPython, File fEnsemble, String strSeqFile, String strQModCheck_mainens, String strQModCheck, String strModCheckPotDat, String strTmJury) throws Exception {
      File fQmo = new File(strTmpPath + strUUID + ".ensemble.qmo");
      RunQModCheck(strQModCheck_mainens, strQModCheck, strModCheckPotDat, strTmpPath, fEnsemble.getCanonicalPath(), fQmo.getCanonicalPath()); //outputs .qmo
      RunTmJury(strTmJury, strTmpPath, strUUID, strSeqFile, fQmo.getCanonicalPath(), fEnsemble.getCanonicalPath(), "10");

      System.out.println("    Building Final Model    ");
      String strCurCmd = strPython + " " + strTmpPath + strUUID + ".jurymod.py";
      System.out.println("    "+strCurCmd);
      ExternalProcess epModPy = new ExternalProcess(strCurCmd.split(" "));
      epModPy.setEnv("PYTHONPATH", strModLib+":"+strModArch);
      epModPy.setEnv("LD_LIBRARY_PATH", strModArch);
      epModPy.setDir(new File(strTmpPath));

      int nRes = epModPy.call();
      Thread.sleep(100);
      if (nRes != 0) {
          //job.UpdateStatus(StatusClass.Error, StatusCode.JobRunning, "psipass2 failed", strError);
          throw new Exception(epModPy.getError());
      }

      //clear up some modeller files
      File f = new File(strTmpPath + strUUID+".ini");
      f.delete();
      f = new File(strTmpPath + strUUID+".rsr");
      f.delete();
      f = new File(strTmpPath + strUUID+".sch");
      f.delete();
      f = new File(strTmpPath + strUUID+".ini");
      f.delete();
      f = new File(strTmpPath + strUUID+".D00000001");
      f.delete();
      f = new File(strTmpPath + strUUID+".V99990001");
      f.delete();
  }



  protected static File catModels (String strTmpPath, String strUUID) throws Exception {
      Pattern pattern = Pattern.compile(strUUID+".+\\.pdb$");

      System.out.println("    Catting PDB files matched to: "+  pattern.pattern());

      File dir = new File(strTmpPath);
      File[] files = dir.listFiles();

      ByteArrayOutputStream baos=new ByteArrayOutputStream();
      for (int i = 0; i < files.length; i++) {
          String fileName = files[i].getName();
          Matcher matcher = pattern.matcher(fileName);
          if (matcher.find())
          //if (fileName.endsWith("pdb"))
          {
              System.out.println("    Catting file: " +  files[i].getName());
              File file = new File(files[i].getPath());
              FileInputStream fis=new FileInputStream(file);
              while(true){
                  int available=fis.available();
                  if(available==0){
                      //System.out.println("0 available");
                    break;
                  }
                  byte[] temp=new byte[available];
                  int read=fis.read(temp);

                  baos.write(temp);
                  if(read==-1){
                    break;
                  }
              }
          }

      }
      File fEnsemble = new File(strTmpPath + strUUID + ".ensemble.pdb");
      Utils.strToFile(fEnsemble.getCanonicalPath(), baos.toString());
      return(fEnsemble);
  }


  protected static String readSequence(String strFile) throws Exception {
    String fastaContents = Utils.fileToStr(strFile);
    fastaContents = fastaContents.substring(fastaContents.indexOf('\n')+1);
    fastaContents = fastaContents.replaceAll("(\\r|\\n)","");
    return(fastaContents);
  }

protected static void tidyGenModels(File fPResults, String strTmpPath, String rootUUID) throws Exception {
        String strRes = Utils.fileToStr(fPResults.getCanonicalPath());
        String[] pResults = strRes.split("\\n");

        int goodHitCount = 0;
        Map<String,String> map = new HashMap<String,String>();
	      for (int i = 0; i < pResults.length; i++)
        {
            if (pResults[i].startsWith("CERT") || pResults[i].startsWith("HIGH") || pResults[i].startsWith("MEDIUM"))
            {
                String[] toks = pResults[i].split("\\s+");
                String pdbID = toks[9];
                map.put(pdbID,pdbID);
                goodHitCount++;
            }
        }

        if (goodHitCount > 0)
        {
            File dir = new File(strTmpPath);
            String[] files = dir.list();
            Pattern pattern = Pattern.compile(rootUUID + "_(.{6}).model.pdb");
            for (int i = 0; i < files.length; i++)
            {
                String file = files[i];
                //System.out.println(file);
                Matcher matcher = pattern.matcher(file);
                if (matcher.find())
                {
                    String matchedPdbId = matcher.group(1);
                    if (!map.containsKey(matchedPdbId))
                    {
                        System.out.println("Deleting : " + strTmpPath + file);
                        File f = new File(strTmpPath + file);
                        f.delete();
                    }
                }
            }
        }
    }

  protected static void runModeller(String strUUID, String strType, String strTmp, String strModellerBin) throws Exception {
    //System.out.println("Pattern " + strUUID+strExt+"."+strType+"_d+.mod.py");
    Pattern pattern = Pattern.compile(strUUID+"\\."+strType+"_(\\d+)\\.mod\\.py");
    System.out.println("    Running MODELLER on PSI-BLAST targets");
    File dir = new File(strTmp);
    String[] files = dir.list();
    if (files == null)
    {
        throw new Exception("Could not list tmp directory");
    }
    else
    {
      //for(int i=0; i<files.length; i++)
      for(int i=0; i<files.length; i++)
      {
        String file = files[i];
        //System.out.println(file);
        Matcher matcher = pattern.matcher(file);
        if (matcher.find())
        {
          String number = matcher.group(1);

          //System.out.println(file);
          String strCurCmd = strModellerBin + " " + strTmp + file;
          System.out.println("    "+strCurCmd);
          ExternalProcess epModPy = new ExternalProcess(strCurCmd.split(" "));
          epModPy.setDir(new File(strTmp));
          int nRes = epModPy.call();
          Thread.sleep(100);

          //clear up tmp modeller files
          File f = new File(strTmp+strUUID+".npj.psiblast_"+number+".mod.log");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".ini");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".rsr");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".sch");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".ini");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".D00000001");
          f.delete();
          f = new File(strTmp+strUUID+"_"+number+".V99990001");
          f.delete();
          //run Modeller on these files
        }
      }
    }
  }

protected static void outputPsiblastModellerInput(String strRes, Double eValueThreshold,
                                             String strTmpPath, String strPDB, String strPDBAA,
                                             String UUID, Double percentIDThreshold) throws Exception {
      List<SeqSimilaritySearchResult> results = parseBlast(strRes);
      //System.out.println("PUSHED BLAST TO LIST\n");
      int alignCount = 0;
      //
      HashMap hashPdbSeqs = new HashMap();
      hashPdbSeqs = readPdbAa(strPDBAA);
      //loop through the list output the alignment to the PIR file then run modeller
      for (SeqSimilaritySearchResult i : results) {
          List<SeqSimilaritySearchHit> hits = i.getHits();

          for (SeqSimilaritySearchHit j : hits) {

              alignCount++;

              if (j.getEValue() < eValueThreshold) {
                  double pcntID = 0.0;
                  String strAlnQ = "", strAlnS = "";
                  // Need to stitch the alignment together from the HSPs
                  // and pad out to cover the target
                  String subjct = j.getSubjectID();
                  String strSubjctSeq = "";

                  if(hashPdbSeqs.containsKey(subjct))
                  {
                      strSubjctSeq = hashPdbSeqs.get(subjct).toString();
                  }

                  List<String> mergedHSPAlignment = stitchBLASTAln(j, strSubjctSeq);

                  strAlnQ = mergedHSPAlignment.get(0);
                  strAlnS = mergedHSPAlignment.get(1);

                  pcntID = Utils.percentID(strAlnQ, strAlnS);

                  // System.out.println("Query"+subjct+"\n");
                  // System.out.println("Query"+strAlnQ+"\n");
                  // System.out.println("Subject"+strAlnS+"\n");
                  if(pcntID >= percentIDThreshold)
                  {
                      printModAlign(subjct.toLowerCase(),strAlnQ,strAlnS,alignCount, strTmpPath,"psiblast",UUID);
                      printModScript(subjct.toLowerCase(),strTmpPath, alignCount, "psiblast", strPDB, UUID);
                  }
              }
          }
      }
  }

  protected static void printModScript(String subjct, String strTmpPath, Integer alignCount, String strType, String strPDB, String uuid) throws IOException {
      String chain = subjct.substring(4, 5);
      String pdbid = subjct.substring(0, 4);
      String strAlignPath = strTmpPath+uuid + "."+strType+"_"+alignCount+".pir";
      String strModScript = "from modeller import *\n";
      strModScript += "from modeller.automodel import *\n\n";
      strModScript += "env = environ()\n";
      strModScript += "env.io.atom_files_directory = [\'"+strPDB+"\',\'"+strTmpPath+"\']\n";
      strModScript += "a = automodel(env, alnfile=\'"+strAlignPath+"\',\n";
      strModScript += "              knowns=(\'"+pdbid+"\', ),\n";
      strModScript += "              sequence=\'"+uuid+"_"+alignCount+"\')\n";
      strModScript += "a.starting_model = 1\n";
      strModScript += "a.ending_model = 1\n";
      strModScript += "a.make()\n";
      strModScript += "ok_models = filter(lambda x: x[\'failure\'] is None, a.outputs)\n";

      strModScript += "key = \'molpdf\'\n";
      strModScript += "ok_models.sort(lambda a,b: cmp(a[key], b[key]))\n";
      strModScript += "m = ok_models[0]\n";
      strModScript += "print \"Top model:%s\" % m['name']\n";

      //File fPDBtarget = new File(strPDB+"pdb"+pdbid+".ent");
      //File fPDBDestination = new File(strTmp+"pdb"+pdbid+".ent");
      //copyFile(fPDBtarget, fPDBDestination);

      String modpyPath = strTmpPath+uuid + "."+strType+"_"+alignCount+".mod.py";
      //System.out.println("Printing PSIBLAST mod.py: "+modpyPath);
      File fModpy = new File(modpyPath);
      Utils.strToFile(fModpy.getCanonicalPath(), strModScript);
    }

  protected static void printModAlign(String subjct, String strAlnQ, String strAlnS, Integer alignCount, String strTmpPath, String strType, String uuid) throws IOException {
      String chain = subjct.substring(4, 5);
      String pdbid = subjct.substring(0, 4);
      String pirAlignment = ">P1;"+uuid+"_"+alignCount+"\n";
      pirAlignment += "sequence:::::::::\n";
      pirAlignment += strAlnQ+"\n";
      pirAlignment += "*\n";
      pirAlignment += ">P1;"+pdbid+"\n";
      pirAlignment += "StructureX:"+pdbid+":FIRST:"+chain.toUpperCase()+":LAST:"+chain.toUpperCase()+"::::\n";
      pirAlignment += strAlnS+"\n";
      pirAlignment += "*\n";
      String pirPath = strTmpPath+uuid + "."+strType+"_"+alignCount+".pir";
      //System.out.println("Printing PSIBLAST PIR: "+pirPath);
      File fPir = new File(pirPath);
      Utils.strToFile(fPir.getCanonicalPath(), pirAlignment);
  }

  public static List<SeqSimilaritySearchResult> parseBlast(String strRes) throws IOException, UnsupportedEncodingException, SAXException {

      // parse blast results....
      BlastLikeSAXParser parser = new BlastLikeSAXParser();
      parser.setModeLazy();

      //make the SAX event adapter that will pass events to a Handler.
      SeqSimilarityAdapter adapter = new SeqSimilarityAdapter();

      //set the parsers SAX event adapter
      parser.setContentHandler(adapter);

      //The list to hold the SeqSimilaritySearchResults
      List<SeqSimilaritySearchResult> results = new ArrayList();

      //create the SearchContentHandler that will build SeqSimilaritySearchResults
      //in the results List
      SearchContentHandler builder = new BlastLikeSearchBuilder(results, new DummySequenceDB("queries"), new DummySequenceDBInstallation());

      //register builder with adapter
      adapter.setSearchContentHandler(builder);

      //parse the file, after this the result List will be populated with
      //SeqSimilaritySearchResults
      ByteArrayInputStream bis = new ByteArrayInputStream(strRes.getBytes("UTF-8"));
      parser.parse(new InputSource(bis));

      return results;
  }


  public static HashMap readPdbAa(String strPDBAA) throws IOException {
      //System.out.println("READ START!\n");
      HashMap hashPdbSeqs = new HashMap();

      FileInputStream fin = new FileInputStream(strPDBAA);
      FileChannel fc = fin.getChannel();
      BufferedReader input = new BufferedReader(new InputStreamReader(fin));
      String line = null;
      String pdbId = null;
      //System.out.println("BEFORE WHILE!\n");
      List<String> lines = new ArrayList<String>();
      while (( line = input.readLine()) != null){
          lines.add(line);
      }

      for (String entry : lines) {
            if(entry.startsWith("\n"))
            {
                continue;
            }

            entry.replace("\n", "");

            if(entry.startsWith(">"))
            {
               pdbId = entry.substring(1,7);
               //System.out.println("MATCHING" + pdbId);
            }
            else
            {
                if(hashPdbSeqs.containsKey(pdbId))
                {
                  String strTmpLine = "";
                  strTmpLine = hashPdbSeqs.get(pdbId).toString();
                  strTmpLine += entry;
                  hashPdbSeqs.put(pdbId, strTmpLine);
                }
                else
                {
                  hashPdbSeqs.put(pdbId, entry);
                }
            }

       }

      return hashPdbSeqs;
  }

  public static List<String> stitchBLASTAln(SeqSimilaritySearchHit hits, String strSubjctSeq) throws IOException {

      List<SimpleSeqSimilaritySearchSubHit> subhits = hits.getSubHits();
      StringBuffer alignedQuerySeq = new StringBuffer(), alignedSubjectSeq = new StringBuffer();

      // We need to make sure they are in order wrt the query
      // sequence

      ArrayList<SeqSimilaritySearchSubHit> shal = new ArrayList<SeqSimilaritySearchSubHit>(subhits);

      Collections.sort(shal, new subHitComparator());

      int current_query_pos = 0;

      for (SeqSimilaritySearchSubHit subhit : subhits) {
          // Add requisite gaps to ensure we have the
          // entire query sequence in the alignment

          for (int i = current_query_pos; i < subhit.getQueryStart() - 1; i++) {
              alignedQuerySeq.append(seq.charAt(i));
              alignedSubjectSeq.append("-");
          }

          Alignment al = subhit.getAlignment();

          List l = al.getLabels();
          String querySubSeq = al.symbolListForLabel(l.get(0)).seqString();
          String subjectSubSeq = al.symbolListForLabel(l.get(1)).seqString();

          for (int i = 0; i < querySubSeq.length(); i++) {
              alignedQuerySeq.append(querySubSeq.charAt(i));
              alignedSubjectSeq.append(subjectSubSeq.charAt(i));
          }
          current_query_pos = subhit.getQueryEnd();
      }

      // Finally add end padding

      for (int i = current_query_pos; i < seq.length(); i++) {
          alignedQuerySeq.append(seq.charAt(i));
          alignedSubjectSeq.append("-");
      }
      ArrayList<String> al = new ArrayList<String>(2);

      //At this point we've padded the Subject/hit with - at the front and end if the Query alignment didn't start at 1 or at the end of the m_seq

      //Now we pad the Query with a bit of seq for anything missing from the Subject/hit.
      current_query_pos = 0;
      int subjectLength = strSubjctSeq.length();
      for (SeqSimilaritySearchSubHit subhit : subhits) {
          if(subhit.getSubjectStart() > 1) {
              //take substring of initial bit of subject append to alignedSubectSeq
              String strInitialSubjct = null;
              try{
                 strInitialSubjct = strSubjctSeq.substring(0, subhit.getSubjectStart()-1);

                  alignedSubjectSeq.insert(0,strInitialSubjct);

                  //generate string of - the right lenfth and append to alignedQuerySeq
                  String strInitialQuery = "";
                  for(int i = current_query_pos; i < subhit.getSubjectStart()-1; i++)
                  {
                      strInitialQuery+="-";
                  }

                  alignedQuerySeq.insert(0,strInitialQuery);

              } catch(Exception ex)
              {
                  strInitialSubjct = "";
              }


          }
          if(subhit.getSubjectEnd() < subjectLength) {
              String strTailingSubject = strSubjctSeq.substring(subhit.getSubjectEnd(), subjectLength);
              alignedSubjectSeq.append(strTailingSubject);

              String strTailingQuery = "";
              for(int i = subhit.getSubjectEnd(); i < subjectLength-1; i++)
              {
                  strTailingQuery+="-";
              }
              alignedQuerySeq.append(strTailingQuery);
          }

      }


      al.add(alignedQuerySeq.toString());
      al.add(alignedSubjectSeq.toString());

      return al;
  }


      private static class subHitComparator implements Comparator<SeqSimilaritySearchSubHit> {

          public subHitComparator() {
          }

          public int compare(SeqSimilaritySearchSubHit A, SeqSimilaritySearchSubHit B) {
              if (A.getQueryStart() < B.getQueryStart()) {
                  return 1;
              } else if (A.getQueryStart() == B.getQueryStart()) {
                  return 0;
              }
              return -1;
          }
      }

}
