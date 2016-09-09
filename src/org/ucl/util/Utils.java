/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.ucl.util;

import java.io.*;
import java.util.HashMap;
import java.util.Vector;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SequenceAlignment;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Alignment;

/**
 *
 * @author ucacmis
 *
 * Just a concatenation of static methods for various purposes
 *
 */
public class Utils {

    public static HashMap<String, StringBuffer> parseAln(String strAln) {

        HashMap<String, StringBuffer>hm_out=new HashMap<String, StringBuffer>();

        String [] lines=strAln.split("\n");

        for(int i=0;i<lines.length;i++){
            if(lines[i].matches("(^\\d\\w{3,5}\\s+.*)|(^target.*)")){
                String [] toks=lines[i].split("\\s+");

                StringBuffer sb=hm_out.get(toks[0]);
                if(sb==null){
                    sb=new StringBuffer();
                    hm_out.put(toks[0], sb);
                }
                sb.append(toks[1]);
            }
        }
        return hm_out;
    }

    public static double percentID(String seqA, String seqB){
        double ret=0.0;
        int Alen=0, Blen=0, match=0;

        for(int i=0;i<seqA.length();i++){
            char charA=seqA.charAt(i);
            char charB=seqB.charAt(i);

            if(charA!='-')
                ++Alen;
            if(charB!='-')
                ++Blen;

            if((charA==charB) && (charA!='-'))
                    ++match;
        }

    // We assume that sequence A is the query

    if(Alen!=0)
        ret=100.0*(double)match/(double)Alen;
    return ret;
    }
    public static double coverage(String seqA, String seqB, int [] cov_arr){
        // Gives us the % coverage of sequence A by sequence B

        double ret=0.0;
        double len=0.0;
        int pos=-1;

        for(int i=0;i<seqA.length();i++){
            char cA=seqA.charAt(i), cB=seqB.charAt(i);
            if(cA!='-'){
                len+=1;
                ++pos;
                if(pos>=cov_arr.length) break;
                if(cB!='-'){
                    ret+=1;
                    ++cov_arr[pos];
                }
            }
        }
    return (ret/len);
    }

    public static String convertPsiPred(String ssPred, String ssConf){
        return "not yet";
    }

    public static int posinseq(String str, int offset){
        int ret=-1, i=0;

        for(i=0;i<str.length();i++){
             char c=str.charAt(i);
             if(c!='-') ++ret;
             if(ret==offset) break;
        }

        // Does this absurd hack work???
        if(i==str.length()) --i;
        return i;
    }

    public static String convertSeqBack(String incoming, String tmplt){
        // Converts the aligned sequence in incoming to the one
        // in tmplt...leaves gaps standing...tmplt must *not* be gapped

        StringBuilder sbl=new StringBuilder(incoming.length());

        int pos=-1;

        for(int i=0;i<incoming.length();i++){
            char c=incoming.charAt(i);
            if(c!='-'){
                char t=tmplt.charAt(++pos);
                c=t;
            }
            sbl.append(c);
        }
    return sbl.toString();
    }

    public static double totalCoverage(int [] cov_arr){
        int tot=0;
        for(int i=0;i<cov_arr.length;i++)
            if(cov_arr[i]>0)
                ++tot;
    return ((double)tot/(double)cov_arr.length);
    }
    public static String deflate(String seq){
        StringBuffer sb=new StringBuffer();
        for(int i=0;i<seq.length();i++){
           char c=seq.charAt(i);
           if(c!='-')
               sb.append(c);
        }
    return sb.toString();
    }
    public static void addSeqToPIR(String strQuery, StringWriter strAlign) {
        int i = 0;
        int nLen = strQuery.length();

        for (i = 0; i * 70 < nLen; i++) {
            String strLine = new String();
            int nEnd = (i + 1) * 70;
            if (nEnd > nLen) {
                nEnd = nLen;
            }

            for (int j = (i * 70); j < nEnd; j++) {
                strLine += strQuery.charAt(j);
            }
            strAlign.write("\n" + strLine);
        }

        strAlign.write("*\n");
    }
    public static String [] reAlign(String strMask, String strTarget) throws Exception {
        Alignment res = null;

        FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");

        SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File("src/BLOSUM62"));

        Sequence query = ProteinTools.createProteinSequence(strMask, "Query");
        Sequence target = ProteinTools.createProteinSequence(strTarget, "Target");

        SequenceAlignment aligner = new NeedlemanWunsch(
                (short)0, // match
                (short)0, // replace
                (short)12, // insert
                (short)12, // delete
                (short)2, // gap extend
                matrix);

        res = aligner.getAlignment(query, target);

        //System.out.println("Smith-Waterman realignment:\n" + aligner.getAlignmentString());

        String seqQuery = res.symbolListForLabel("Query").seqString();
        String seqSubject = res.symbolListForLabel("Target").seqString();

        String [] out = new String [2];

        out[0]=seqQuery.replace("~", "-");
        out[1]=seqSubject.replace("~", "-");

        return out;

    }

    public static double pcntID(String strA, String strB){
        if(strA.length()!=strB.length()) return -1.0;

        int lenA=0, lenB=0;
        int same=0;

        for(int i=0;i<strA.length();i++){
            char q=strA.charAt(i);
            char s=strB.charAt(i);

            if(q!='-')
                ++lenA;
            if(s!='-')
                ++lenB;

            if(q!='-' && (q==s))
                ++same;
        }

        int len=lenA<lenB?lenA:lenB;

        if(len==0) len=1;

        return 100.0*(double)same/(double) len;
    }

    public static String clipGappedAlignment(String strMask, String strTarget) {
        String strRes = new String();
        if (strMask.length() != strTarget.length() || strMask.isEmpty()) {
            return null;
        }
        int nGapLen = 3;
        // replaces all masked elements with gap size > nGapLen with gap symbols '-'

        // lazy, scan gaps, then mask them
        int nIn = 0;
        int nOut = 0;

        int i = 0;
        boolean bGap = false; // false for regular, true for gap.
        Vector<Integer> vGaps = new Vector<Integer>();

        // make a fast gapmap that satisfies gaplen,
        // NOTE: not DP, just simple state machine

        if (strMask.charAt(0) == '-') {
            nIn = 0;
            bGap = true;
        }
        for (i = 0; i < strMask.length(); i++) {
            if (bGap) {
                if (strMask.charAt(i) != '-') {
                    bGap = false;
                    if ((i - nIn) > nGapLen) {
                        vGaps.add(nIn);
                        vGaps.add(i);
                    }
                }
            } else {
                if (strMask.charAt(i) == '-') {
                    bGap = true;
                    nIn = i;
                }
            }
        }
        // short circuit if no gaps to mask
        if (vGaps.isEmpty()) {
            return strTarget;
        }

        int nGapCount = vGaps.size() / 2;
        int j = -1;
        nIn = vGaps.elementAt(0);
        nOut = vGaps.elementAt(1);

        for (i = 0; i < strTarget.length(); i++) {
            if ((i > nOut) && (j < nGapCount - 1)) {
                j++;
                nIn = vGaps.elementAt(j * 2);
                nOut = vGaps.elementAt((j * 2) + 1);
            }

            if ((i >= nIn) && (i < nOut)) {
                strRes = strRes + "-"; // mask it
            } else {
                strRes = strRes + strTarget.charAt(i);
            }
        }

        return strRes;
    }

    public static String getModellerAtomSequence(String pdbID, String temp_path,
                                                 String pdb_path, String modeller_bin)throws Exception{

        // Find the requisite chain

        String chain=pdbID.substring(4,5);

        String strModPy="from modeller import *\n";
        strModPy+="env = environ()\n";
        strModPy += "env.io.atom_files_directory = \'"+pdb_path+"\'\n" ;
        strModPy+="code = \'"+ pdbID.substring(0,4) +"\'\n";
        strModPy+="mdl = model(env, file=code,model_segment=(\'FIRST:" + chain + "\', \'LAST:" + chain +"\'))\n";
        strModPy+="aln = alignment(env)\n";
        strModPy+="aln.append_model(mdl, align_codes=code)\n";
        strModPy+="aln.write(file=\'" +temp_path +"\'+code+\'.seq\')\n";
        File fModPY= new File(temp_path+ "_forseq.py");
        strToFile(fModPY.getCanonicalPath(), strModPy);
        String strCurCmd = modeller_bin+ " " + fModPY.getCanonicalPath();
        ExternalProcess epModPy = new ExternalProcess(strCurCmd.split(" "));
        System.out.println("RUNNING MODELLER WITH: " + strCurCmd);

        try
        {
            int nRes = epModPy.call();
            Thread.sleep(100); // let output catch up =(
        }
        catch (Exception ex)
        {
            System.out.println("MODELLER ISSUE" + epModPy.getError());
        }
        //modeller generates this .seq file
        System.out.println("MODELLER OUTPUT SEQ: " + pdbID.substring(0,4));
        String seqOutput=temp_path+pdbID.substring(0,4)+".seq";

        String strOutput=fileToStr(seqOutput);

        String [] seqs=strOutput.split("\n");


        boolean gettit=false;

        String strSeq="";

        for(int i=0;i<seqs.length;i++){
            //System.out.println("MODELLER OUTPUT SEQ: " +seqs[i]);
            if(seqs[i].contains("structure")){
                gettit=true;
                continue;
            }
            if(gettit && !seqs[i].isEmpty())
                strSeq+=seqs[i];
        }

        //fModPY.delete();
        //new File(seqOutput).delete();
        return strSeq.replace("*", "");
    }

    public static String addSS(String strSS) throws Exception {
        String[] entries = strSS.split("\n");

        StringBuffer sb=new StringBuffer();

        for (int i = 2; i < entries.length; i++) {
            String[] ssRec = entries[i].trim().split("\\s+");

            char cType=ssRec[2].charAt(0);
            float fScore=0;

            float fCoil = Float.parseFloat(ssRec[3].trim());
            float fHelix = Float.parseFloat(ssRec[4].trim());
            float fStrand = Float.parseFloat(ssRec[5].trim());



            if (cType == 'C') {
                fScore = fCoil;
            } else if (cType == 'H') {
                fScore = fHelix;
            } else if (cType == 'E') {
                fScore = fStrand;
            }

            // davids classifier
            /*#define MIN(x,y) (((x)<(y))?(x):(y))
            #define MAX(x,y) (((x)>(y))?(x):(y))
            main(int argc, char **argv)
            {
            char buf[512];
            float cp, hp, ep, cutoff = 0.8;
            if (argc > 1)
            cutoff = atof(argv[1]);
            while (fgets(buf, 160, stdin))
            {
            if (sscanf(buf+8, "%f%f%f", &cp, &hp, &ep) != 3)
            continue;
            if (2*MAX(MAX(cp, hp), ep)-(cp+hp+ep)+MIN(MIN(cp, hp), ep) > cutoff)
            putchar(buf[7]);
            else if (cp < MIN(ep, hp))
            putchar('?');
            else if (ep > hp)
            putchar('e');
            else if (hp > ep)
            putchar('h');
            else
            putchar('?');
            }
            putchar('\n');
            } */
            float fCutoff = 0.8f;
            if ( ( (2*fScore) -
                    ( (fCoil+fStrand+fHelix)
                    + (Math.min(Math.min(fCoil, fStrand), fHelix) ) ) ) <= fCutoff )
            {
                if (fCoil < Math.min(fStrand, fHelix))
                {
                    cType = '?';
                } else if (fStrand > fHelix)
                {
                    cType = 'e';
                } else if (fHelix > fStrand)
                {
                    cType = 'h';
                } else
                {
                    cType = '?';
                }
            }

            sb.append(cType);
       }
      return sb.toString();
    }

    public static String fileToStr(String strInput) throws IOException {

        // It is implicit that the input is a file so...

        File f=new File(strInput);

        // Ensure we have the capacity to prevent constant reallocation

        StringBuilder sb=new StringBuilder(new Long(f.length()).intValue());

        FileInputStream fIn = new FileInputStream(strInput);
        BufferedReader br = new BufferedReader(new InputStreamReader(fIn));
        String line = null;
        while ((line = br.readLine()) != null) {
            sb.append(line);
            sb.append("\n");
        }

        br.close();

        return sb.toString();
    }

    public static void strToFile(String strPath, String strOutput) throws IOException {
        //System.out.println(strOutfile+"\n");
        FileOutputStream fOut = new FileOutputStream(strPath);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fOut));
        bw.write(strOutput);
        bw.flush();
        bw.close();
    }
}
