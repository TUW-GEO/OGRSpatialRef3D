/*
 * $Id: ogr_php.i 8189 2005-09-02 16:19:23Z kruland $
 *
 * php specific code for ogr bindings.
 */

/*
 * $Log$
 * Revision 1.1  2005/09/02 16:19:23  kruland
 * Major reorganization to accomodate multiple language bindings.
 * Each language binding can define renames and supplemental code without
 * having to have a lot of conditionals in the main interface definition files.
 *
 */

%init %{
  if ( OGRGetDriverCount() == 0 ) {
    OGRRegisterAll();
  }
%}

%include typemaps_php.i
