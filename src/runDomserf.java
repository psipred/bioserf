class runDomserf
{
  static String strJobExt = "bioserf";
  static String seq;

  public static void main( String args[] ) throws Exception
  {
    if( args.length != 1 )
    {
          System.out.println( "usage: java -cp bioserf/src/:bioserf/src/org/ucl/ runDomserf B0R5N0.domserf.bls" );
          System.exit(1);
    }
  String strBlastFile = args[0]; //the input blast file, blast vs CathDomainSeqs.S100.ATOM
  System.out.println(strBlastFile);
  }

}
