/*
 * system.i --
 *
 *	System related utilities for Yorick.  Provides routines:
 *	  exec  - similar to system routine but allow command splitted in
 *	          multiple pieces and allow to get the command result.
 *	  mktemp - gererate temporary file name.
 *	  get_filesize - get size of file.
 *	  remove_log_file - remove Yorick's contents log file.
 *
 * History:
 *	$Id: system.i,v 1.3 1997/04/03 15:40:56 eric Exp $
 *	$Log: system.i,v $
 *	Revision 1.3  1997/04/03 15:40:56  eric
 *	Fixed typo error.
 *
 *	Revision 1.2  1997/04/03 15:39:14  eric
 *	Make use of the new Yorick's command popen() to avoid temporary files.
 *
 *	Revision 1.1  1996/08/05 09:29:13  eric
 *	Initial revision
 *
 */

/*---------------------------------------------------------------------------*/
func remove_log_file(name, interactive=)
/* DOCUMENT ok = remove_log_file(name, interactive=1/0)
     Remove Yorick's log file NAME+"L".  If INTERACTIVE is non-nil and
     non-zero, user is asking for confirmation.  If a log file exists
     and is not removed, 0n is returned; otherwise, 1n is returned.
 */
{
    nameL = name + "L";
    if (open(nameL, "", 1)) {
	if (interactive) {
	    write, format="%s", "Remove existing file "+nameL+" (y or n)? ";
	    answer = strtok(rdline(prompt=string(0)))(1);
	    if (numberof(where(answer == ["y", "Y", "yes", "Yes", "YES"]))
		    == 0) {
		write, "No file written or deleted.";
		return 0n;
	    }
	}
	remove, nameL;
    }
    return 1n;
}
/*---------------------------------------------------------------------------*/
func mktemp(template, dir=)
/* DOCUMENT tmpname = mktemp(template, dir=)
     Search for a file name not yet in use.  If TEMPLATE is given, it is
     used as a template to generate temporary file name.  If TEMPORARY
     contains the string "%d", it will be substituted by a number; otherwise
     the number will be appended to TEMPLATE.  If TEMPLATE is omitted, it
     defaults to: "yorick.tmp%d".

     Keyword DIR can be used to specify the directory into which
     to create the temporary file.  By default DIR is (in order):
       (1) an empty string if TEMPLATE starts with a "/" or a "./" or
           a "../";
       (2) the contents of environment variable TMPDIR if it exists;
       (3) "/tmp". 
     
   RESTRICTIONS:
     (1) There is no way to be sure that a file with the same name will not
         be created after the moment when the temporary file name is checked
         to not yet exists.
     (2) Do not forget to remove the temporary file when no more needed; use
         Yorick's subroutine "remove" for that.  This must be done even if
	 the temporary file is not used, since mktemp creates an empty
	 file TMPNAME to reduce problems as explained in (1).
     (3) If TEMPLATE contains some "/", the effective directory where will
         be the temporary file may not be DIR.  The effective directory
	 must be writable.  An error is raised if this not the case since
	 mktemp creates an empty file TMPNAME.
     
   SEE ALSO:
     exec, open, remove, get_env.
 */
{
    if (is_void(template)) {
	template = "yorick.tmp%d";
    } else if (numberof(where(strmatch(template, "%d"))) == 0) {
	template += "%d";
    }
    if (is_void(dir)) {
	if (strpart(template, 1:1) != "/" && strpart(template, 1:2) != "./"
		&& strpart(template, 1:3) != "../") {
	    dir = get_env("TMPDIR");
	    if (!dir)
		dir = "/tmp";
	    template = dir + "/" + template;
	}
    } else {
	template = dir + "/" + template;
    }
    i = 0n;
    do {
	name = swrite(format=template, i++);
    } while (open(name, "", 1));
    open, name, "w";	// Create an empty temporary file.
    return name;
}
/*---------------------------------------------------------------------------*/
func exec(cmd, ..)
/* DOCUMENT exec, shell_command, ...
            result = exec(shell_command, ...)
     Concatenates its arguments with spaces inserted between them to form a
     shell command string.  If called as a function, returns command output
     as an array of strings; otherwise, just execute the command.

   SEE ALSO: system, popen.
 */
{
    while (more_args())
	cmd += " " + next_arg();
    if (am_subroutine()) {
	system, cmd;
	return;
    }
    pipe= popen(cmd, 0);
    result= rdline(pipe, 1000);
    while (result(0))
	grow, result, rdline(pipe, 1000);
    return result(1:where(!result)(1)-1);
}
/*---------------------------------------------------------------------------*/
func get_filesize(name)
/* DOCUMENT nbytes = get_filesize(name)
     Returns size of NAME in bytes.

   RESTRICTIONS:
     Makes use of UNIX command "wc" which may not exists on other systems.

   SEE ALSO: exec.
 */
{
    if (!open(name, "", 1))
	error, "no such file or directory \""+name+"\"";
    sz= 0L;
    if (sread(rdline(popen("wc -c " + name, 0)), sz) != 1)
	error, "cannot get size of file \""+name+"\"";
    return sz;
}
/*---------------------------------------------------------------------------*/
