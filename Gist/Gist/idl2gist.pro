PRO idl2gist, filename
;+
; NAME:
;	idl2gist
;
; PURPOSE:
;	Convert IDL color tables to Yorick-Gist palettes.
;
; CALLING SEQUENCE:
;	idl2gist [, inputFileName]
;
; SIDE EFFECTS:
;	The files idl-01.gp, idl-02.gp, ..., are created.
;
; RESTRICTIONS:
;	Works from the file: $IDL_DIR/resource/colors/colors1.tbl or the file specified
;	with the FILE keyword.
;
; PROCEDURE:
;	The file "colors1.tbl" or the user-supplied file is read.  If
;       the currently selected device doesn't have 256 colors, the color
;	data is interpolated from 256 colors to the number of colors
;	available.
;
;	The colors loaded into the display are saved in the common
;	block COLORS, as both the current and original color vectors.
;
;	Interpolation:  If the current device has less than 256 colors,
;	the color table data is interpolated to cover the number of
;	colors in the device.
;
; MODIFICATION HISTORY:
;
;-

    on_ioerror, bad
    on_error, 2                 ;Return to caller if error
    get_lun, input

    
    IF (n_elements(filename) GT 0) THEN BEGIN
        filename = string(file)
    END ELSE BEGIN
        filename = filepath('colors1.tbl', subdir=['resource', 'colors'])
    END

    openr,input, filename, /block
    
    ntables = 0b
    readu, input, ntables

    ;; Read names
    names = bytarr(32, ntables)
    point_lun, input, ntables * 768L + 1L ;Read table names
    readu, input, names
    names = strtrim(names, 2)

    ;; There is a maximum of 240 colors in Gist/Yorick
    ncolors= 240L

    FOR table_number=0L, ntables-1L DO BEGIN

        get_lun, output
        palette_name = string(format='("idl-",i2.2,".gp")', table_number)
        idl_name = names(table_number)
        openw, output, palette_name
        
        printf, output, "# Gist " + palette_name + " palette"
        printf, output, "# from IDL "+ idl_name + " palette"
        printf, output
        printf, output, string(format='("ncolors=", i)', ncolors)
        printf, output
        printf, output, "# ntsc gray scale looks slightly better than straight intensity"
        printf, output, "ntsc= 1"
        printf, output
        printf, output, "#  r   g   b"
        
        aa=assoc(input, bytarr(256),1) ;Read 256 long ints
        r = aa(table_number*3L)
        g = aa(table_number*3L+1L)
        b = aa(table_number*3L+2L)
    
        IF (ncolors NE 256L) THEN  BEGIN
            ;; Interpolate
            u = replicate(1, 256) # ((255./ (ncolors-1)) * findgen(ncolors))
            v = findgen(256) # replicate(1, ncolors)
            k = 1. - (abs(u - v) < 1.)
            r = long(r # k + .5) < 255L
            g = long(g # k + .5) < 255L
            b = long(b # k + .5) < 255L
        END
    
        FOR i=0L, ncolors-1 DO printf, output, format='(i4,i4,i4)', r(i), g(i), b(i)
        
        ;; close, output
        free_lun, output        ; also close file?
    end

    free_lun, input             ; also close file?
    return
    
bad:
    message,'Error reading file: ' + filename + ', ' + !err_string
    
end
