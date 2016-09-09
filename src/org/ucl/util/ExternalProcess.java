/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.ucl.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author sward
 */
public class ExternalProcess implements Callable {
    protected String m_strCommandLine[];
    protected StringWriter m_results, m_error;
    protected ProcessBuilder m_pb;
    protected int m_nResCode;
    protected Vector<File> m_fInputs;
    protected String m_strInput;
    
    // To allow easy classification of inputs for remote invokeations
    public void addInput(File fInput) {
        m_fInputs.add(fInput);
    }
    
    public void addStrInput(String strInput){
        m_strInput=strInput;
    }
            
    
    public int getResCode() {
        return m_nResCode;
    }
    
    public String getOutput() {
        if (m_results != null) {
            return m_results.toString();
        } else {
            return null;
        }
    }
    
    public String getError(){
        if(m_error!=null){
            return m_error.toString();
        }
        else{
            return null;
        }
    }

    public void setEnv(String strKey, String strVal)
    {
        Map<String, String> env = m_pb.environment();
        env.put(strKey, strVal);
    }
    
    public void setDir(File fDir)
    {
        m_pb.directory(fDir);
    }
    
    protected class Gobbler extends Thread {
        private InputStream m_is;
        private Writer m_os;
        
        // Use null writer to simply eat the stream
        public Gobbler(InputStream is, Writer os) {
            m_is = is;
            m_os = os;
        }
        
        public void cleanup() {
            try {
                m_is.close();
            } catch (IOException ex) {
                // ignore it, since cleanup is occuring anyways...
            }
        }
        
        // Responsible for eating output from the running process
        // On a cancel, will clean up and pass along the cancel...
        public void run() {
            try {

                String line = null;
                BufferedReader br = new BufferedReader(new InputStreamReader(m_is));
                // eat the stream...
                while ((line = br.readLine()) != null) {
                    if (m_os != null) m_os.write(line + "\n");
                }
            } catch (IOException ex) {
                Logger.getLogger(ExternalProcess.class.getName()).log(Level.SEVERE, null, ex);
            } 
            cleanup();
        }
    }
    
    public ExternalProcess( String strCommandLine[]) {
        m_strCommandLine = strCommandLine;
        m_pb = new ProcessBuilder(m_strCommandLine);
        m_nResCode = -1;
        m_fInputs = new Vector<File>();
        m_strInput=null;
    }
    
    public Integer call() {
        m_nResCode = -1;
        Process p = null;
        try {
            m_results = new StringWriter();
            m_error = new StringWriter();
            
            p = m_pb.start();
            
            if(m_strInput!=null){
                PrintWriter pw=new PrintWriter(p.getOutputStream());
                String [] toks=m_strInput.split("\n");
                for(int i=0;i<toks.length;i++)
                    pw.println(toks[i]);
                pw.close();
            }
                        
            p.getOutputStream().close();
            Gobbler gIn = new Gobbler(p.getInputStream(), m_results);
            Gobbler gErr = new Gobbler(p.getErrorStream(), m_error);
            // TODO: cleaner thread shutdown...
            gIn.start();
            gErr.start();
            
            Thread.yield();
            m_nResCode = p.waitFor();

            Thread.sleep(10);
            
        } catch (IOException ex) {
            Logger.getLogger(ExternalProcess.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(ExternalProcess.class.getName()).log(Level.SEVERE, null, ex);
            p.destroy();
            // killing process will also reach the gobblers...
            // could add explicit kills for them, however.
        }
        return m_nResCode;
    }
}
