"""
  BioPython GUI
  JLG 31/12/2021
 
"""
import sys, os, time
import subprocess
import tkinter as tk
from io import StringIO
from ttkthemes  import ThemedTk
import tkinter as tk
from tkinter    import ttk
from tkinter    import messagebox
from tkinter    import scrolledtext
from tkinter.filedialog import askopenfilename, asksaveasfilename , askdirectory
from Bio import __version__ as bioversion
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import DialignCommandline
from Bio.Align.Applications import ProbconsCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio import Entrez
from Bio import ExPASy
from Bio import SwissProt
from Bio import TogoWS
from Bio.ExPASy import ScanProsite
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqUtils import GC
from Bio.SeqUtils import molecular_weight
import webbrowser

# User Params
Entrez.email = "******@******.**"
PATH_FIREFOX ='C://Program Files//Mozilla Firefox//firefox.exe'

# User Functions
def center_window(w=1420, h=640):
    ws = window.winfo_screenwidth()
    hs = window.winfo_screenheight()
    x = (ws/2) - (w/2)    
    y = (hs/2) - (h/2) - 40
    window.geometry('%dx%d+%d+%d' % (w, h, x, y))
    
def About():
    tk.messagebox.showinfo(title="About ...",
                           message="Python Version:\n" + sys.version +
                                "\nTkinter Version: "  + str(tk.TkVersion) +
                                "\nBiopython Version: "  + bioversion
                           )
    
def Close():
    res = messagebox.askquestion('Exit', 'Do you want to exit?')
    if res == 'yes':
        window.destroy()
        
class CrearToolTip(object):
    def __init__(self,elemento,texto='Info del objeto'):
        self.espera = 500
        self.largo = 180
        self.objeto = elemento
        self.texto = texto
        self.objeto.bind("<Enter>",self.entrar)
        self.objeto.bind("<Leave>",self.salir)
        self.objeto.bind("<ButtonPress>",self.salir)
        self.id = None
        self.tw = None
    def entrar(self,event=None):
        self.asignar()
    def salir(self,event=None):
        self.liberar()
        self.ocultar_tip()
    def  asignar(self):
        self.liberar()
        self.id = self.objeto.after(self.espera,self.mostrar_tip)
    def  liberar(self):
        id = self.id
        self.id = None
        if id:
            self.objeto.after_cancel(id)
    def  mostrar_tip(self,event=None ):
        x = y = 0
        x,y,cx,cy = self.objeto.bbox("insert")
        x += self.objeto.winfo_rootx() + 25
        y += self.objeto.winfo_rooty() + 20
        self.tw = tk.Toplevel(self.objeto)
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d"%(x,y))
        label = tk.Label(self.tw,text=self.texto,justify="left",
                       background = "white",relief="solid",borderwidth=1,
                       wraplength = self.largo)
        label.pack(ipadx=1)
    def ocultar_tip(self):
        tw = self.tw
        self.tw = None
        if tw:
            tw.destroy()
            
def rClicker(e):
    #print (e.widget.message)
    try:
        def rClick_Copy(e, apnd=0):
            e.widget.event_generate('<Control-c>')
        def rClick_Cut(e):
            e.widget.event_generate('<Control-x>')
        def rClick_Paste(e):
            e.widget.event_generate('<Control-v>')
        def rClick_Select_All(e):
            e.widget.focus_force()
            e.widget.tag_add(tk.SEL, "1.0", tk.END)
        def rClick_Clear_All(e):
            e.widget.focus_force()
            e.widget.tag_add(tk.SEL, "1.0", tk.END)
            e.widget.event_generate('<Control-x>')
        def rClick_Clipboard(e):
            window.clipboard_append(e.widget.get(1.0, tk.END ))
        def rClick_To_File(e):           
            filepath = asksaveasfilename(
                 defaultextension="",initialfile ="seq000.fasta",
                 filetypes=[("Seq Files", "*.fasta;*.gb;*.genbank"), ("All Files", "*.*")],
                 )
            if not filepath:
                return
            f = open(filepath, 'w')
            f.write(e.widget.get(1.0, tk.END ))
            f.close           
        def rClick_AliView_seq(e):
            FileView=seq_inuse_entry.get()                    
            subprocess.run(["ExtApp\Aliview.exe", FileView])
        def rClick_AliView_aln(e):              
            FileView=os.path.splitext(seq_inuse_entry.get())[0]+'.aln'
            subprocess.run(["ExtApp\Aliview.exe", FileView])
        def rClick_ClustalX_seq(e):
            FileView=seq_inuse_entry.get()                    
            subprocess.run(["ExtApp\clustalx.exe", FileView])
        def rClick_ClustalX_aln(e):              
            FileView=os.path.splitext(seq_inuse_entry.get())[0]+'.aln'
            subprocess.run(["ExtApp\clustalx.exe", FileView])

        e.widget.focus()
         
        if (e.widget.message=='seq_edit'):
            nclst=[
               ('Cut'  , lambda e=e: rClick_Cut(e)),
               ('Copy' , lambda e=e: rClick_Copy(e)),
               ('Paste', lambda e=e: rClick_Paste(e)),
               ('Select All',        lambda e=e: rClick_Select_All(e)),
               ('Clear All',         lambda e=e: rClick_Clear_All(e)),
               ('Copy to Clipboard', lambda e=e: rClick_Clipboard(e)),
               ('Save to File', lambda e=e: rClick_To_File(e)),
               ('-', lambda e=e: rClick_Clipboard(e)),
               ('Seq with AliView', lambda e=e: rClick_AliView_seq(e)),
               ('Seq with ClustalX', lambda e=e: rClick_ClustalX_seq(e)), 
               ]
        elif (e.widget.message=='aln_edit'):
            nclst=[
               ('Cut'  , lambda e=e: rClick_Cut(e)),
               ('Copy' , lambda e=e: rClick_Copy(e)),
               ('Paste', lambda e=e: rClick_Paste(e)),
               ('Select All',        lambda e=e: rClick_Select_All(e)),
               ('Clear All',         lambda e=e: rClick_Clear_All(e)),
               ('Copy to Clipboard', lambda e=e: rClick_Clipboard(e)),
               ('Save to File', lambda e=e: rClick_To_File(e)),
               ('-', lambda e=e: rClick_Clipboard(e)),
               ('Aln with AliView', lambda e=e: rClick_AliView_aln(e)),
               ('Aln with ClustalX', lambda e=e: rClick_ClustalX_aln(e)),
               ]
        else:
            nclst=[
               ('Cut'  , lambda e=e: rClick_Cut(e)),
               ('Copy' , lambda e=e: rClick_Copy(e)),
               ('Paste', lambda e=e: rClick_Paste(e)),
               ('Select All',        lambda e=e: rClick_Select_All(e)),
               ('Clear All',         lambda e=e: rClick_Clear_All(e)),
               ('Copy to Clipboard', lambda e=e: rClick_Clipboard(e)),
               ('Save to File', lambda e=e: rClick_To_File(e)),
               ]
            
        rmenu = tk.Menu(None, tearoff=0, takefocus=0)
        for (txt, cmd) in nclst:
            rmenu.add_command(label=txt, command=cmd)
        rmenu.tk_popup(e.x_root+40, e.y_root+10,entry="0")
    except tk.TclError:       
        print (' - rClick menu, something wrong')
        pass
    return "break"

def Load_Seq():
    filepath = askopenfilename(filetypes=[("SEQ Files", "*.fasta;*.fas;*.gb;*.gbk;*.genbank"), ("All Files", "*.*")])
    if not filepath:
        return
    f = open(filepath, 'r')
    f2 = f.read()
    seq_edit.delete('1.0', 'end')
    seq_edit.insert('1.0', f2)
    seq_inuse_entry.delete(0, tk.END)
    seq_inuse_entry.insert(0, (filepath))
    f.close()
    
    seq_format=os.path.splitext(filepath)[1]
    #print(seq_format)
    if  ( seq_format in  ['.gbk','.gb'] ) :
        format_seq='genbank'
    elif (seq_format in  ['.fasta','.fas','.fa']):
        format_seq='fasta'
    elif (seq_format in  ['.nexus','.nex','.ne']):
        format_seq='nexus'
    else:
        print('Format not defined ...') 
    label_var.set(format_seq)
    
def Save_Seq():    
    filepath = asksaveasfilename(
                 defaultextension="",initialfile ="seq000.fasta",
                 filetypes=[("Seq Files", "*.fasta;*.gb;*.genbank"), ("All Files", "*.*")],
            )
    if not filepath:
        return
    f = open(filepath, 'w')
    f.write(seq_edit.get(1.0, tk.END))
    f.close
    
def Convert_Seq():
    seq_file=seq_inuse_entry.get()
    format_seq=label_var.get()
    format_dest=Convert_Combo.get()
    conv_edit.delete('1.0', 'end')
    try:
        if  (format_dest=='Fasta'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.fasta'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "fasta")
        elif (format_dest=='Fasta2'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.fasta2'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "fasta-2line")
        elif (format_dest=='Fastq'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.fastq'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "fastq")
        elif (format_dest=='Fastq Sanger'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.fastq'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "fastq-sanger")
        elif (format_dest=='Genbank'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.gb'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "genbank")
        elif (format_dest=='embl'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.embl'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "embl")
        elif (format_dest=='Nexus'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.nexus'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "nexus")
        elif (format_dest=='tab'):
            conv_file=aln_file=os.path.splitext(seq_file)[0]+'.tab'
            count = SeqIO.convert(seq_file, format_seq, conv_file, "tab")  
        else:
            messagebox.showwarning(title="Format Conversion", message='Method not defined ...' )
    except ValueError as v_e:
        messagebox.showwarning(title="Format Conversion", message='Error executing aligment !' + str(v_e)) 
        return False           
    f = open(conv_file, 'r')
    f2 = f.read()   
    conv_edit.insert('1.0', f2)
    messagebox.showwarning(title="Format Conversion", message=' Converted ' + str(count) + ' records' )     
    
def Align_Seq():
    seq_file=seq_inuse_entry.get()
    aln_file=os.path.splitext(seq_inuse_entry.get())[0]+'.aln'
    aln0_file=os.path.splitext(seq_inuse_entry.get())[0]
    aln1_file=os.path.splitext(seq_inuse_entry.get())[0]+'.best.nex'
    tre_file=os.path.splitext(seq_inuse_entry.get())[0]+'.dnd'
    tre0_file=os.path.splitext(seq_inuse_entry.get())[0]+'.best.dnd'
    aln_edit.delete('1.0', 'end')
    dnd_edit.delete('1.0', 'end')
    stdout_edit.delete('1.0', 'end')
    stderr_edit.delete('1.0', 'end')
    label0_var.set(' Align in process ...')
    window.update_idletasks()
    
    if (Align_Combo.get()=='Clustalw'):
        file = ["Align/clustalw2.exe", "-infile="+seq_file, "-TYPE="+SeqType_Combo.get()]
        result=subprocess.run(file, capture_output=True, text=True)
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)            
        tree = Phylo.read(tre_file, "newick")   
        temp_out = StringIO()
        sys.stdout = temp_out
        Phylo.draw_ascii(tree)           
        dnd_edit.insert('1.0', temp_out.getvalue())
        sys.stdout = sys.__stdout__
        Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
        label0_var.set('')
        window.update_idletasks()
        
    
    elif (Align_Combo.get()=='ClustalO'):
        file = ["Align/clustalo.exe", "--infile="+seq_file,"--outfile="+aln_file,"-v","--force","--outfmt=clu"]
        result=subprocess.run(file, capture_output=True, text=True)
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)                  
        label0_var.set('')
        window.update_idletasks()
    
    elif (Align_Combo.get()=='Prank'):
        file = ["Align/prank.exe", "-d="+seq_file,"-o="+aln0_file,"-f=17","-showtree"]
        result=subprocess.run(file, capture_output=True, text=True)
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln1_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)     
        tree = Phylo.read(tre0_file, "newick")   
        temp_out = StringIO()
        sys.stdout = temp_out
        Phylo.draw_ascii(tree)           
        dnd_edit.insert('1.0', temp_out.getvalue())
        sys.stdout = sys.__stdout__
        Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
        label0_var.set('')
        window.update_idletasks()
    
    elif (Align_Combo.get()=='Muscle0'):
        #cmd = "muscle5.exe -align " + seq_file +" -output "+aln_file
        file = ["muscle5.exe", "-align " + seq_file,"-output "+aln_file]     
        proc = subprocess.Popen(file, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        print ("program output:", out)       
        #result=os.popen(cmd).read()    
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)
        f.close()     
        label0_var.set('')
        window.update_idletasks() 
    elif (Align_Combo.get()=='Muscle'):
        try:
            muscle_exe = r"muscle.exe"        
            muscle_cline = MuscleCommandline(input=seq_file, out=aln_file, clw=True, tree2=tre_file )
            print (muscle_cline)
            assert os.path.isfile(muscle_exe), "Muscle executable missing"
            stdout, stderr = muscle_cline()   
            f = open(aln_file, 'r')
            f2 = f.read()        
            aln_edit.insert('1.0', f2)
            f.close()
            tree = Phylo.read(tre_file, "newick")   
            temp_out = StringIO()
            sys.stdout = temp_out
            Phylo.draw_ascii(tree)           
            dnd_edit.insert('1.0', temp_out.getvalue())
            sys.stdout = sys.__stdout__
            Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
            label0_var.set('')
            window.update_idletasks()          
        except BaseException as e:
            messagebox.showwarning(title="Muscle", message=str(e) )
            label0_var.set('')
            window.update_idletasks()
            return False
    elif (Align_Combo.get()=='Dialign'):
        file = ["Align/dialign-tx.exe",  "-D","-l0","-t4","-n4","-g40","mat",seq_file ,aln_file]      
        result=subprocess.run(file, capture_output=True, text=True)
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)                  
        label0_var.set('')
        window.update_idletasks()    
    
    elif (Align_Combo.get()=='Probcons'):
        file = ["Align/probcons.exe", seq_file ,"-clustalw","-v"]      
        result=subprocess.run(file, capture_output=True, text=True)       
        aln_edit.insert('1.0', result.stdout)
        stdout_edit.insert('1.0', result.stderr)
        #stderr_edit.insert('1.0', result.stderr)
        label0_var.set('')
        window.update_idletasks()
        
    elif (Align_Combo.get()=='Kalign'):
        #file = ["Align/kalign.exe", "-i "+seq_file ,"-o "+aln_file, "-b NJ","-f Clu","-c input","-d wu","-s 80","-t 2","-e 3","-m 28.3"]      
        file = ["kalign.exe", seq_file ]         
        result=subprocess.run(file)
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)                  
        label0_var.set('')
        window.update_idletasks()
        
    elif (Align_Combo.get()=='PCMA'):
        file = ["Align/pcma.exe", seq_file ]
        result=subprocess.run(file, capture_output=True, text=True) 
        stdout_edit.insert('1.0', result.stdout)
        stderr_edit.insert('1.0', result.stderr)
        f = open(aln_file, 'r')
        f2 = f.read()         
        aln_edit.insert('1.0', f2)            
        tree = Phylo.read(tre_file, "newick")   
        temp_out = StringIO()
        sys.stdout = temp_out
        Phylo.draw_ascii(tree)           
        dnd_edit.insert('1.0', temp_out.getvalue())
        sys.stdout = sys.__stdout__
        Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
        label0_var.set('')
        window.update_idletasks() 
    
    else:
        messagebox.showwarning(title="Align", message='Method not defined ...' )
        
    label0_var.set('')
    window.update_idletasks()
    
def Entrez_Search():
    entrez_edit.delete('1.0', 'end')
    try:
        handle = Entrez.efetch(db=db_Combo.get(), id=id_entry.get() , rettype=type_entry.get(), retmode=mode_entry.get())
    except BaseException as e:
        messagebox.showwarning(title="Entrez", message=str(e) )
        return False
    mess=handle.read()
    if (mess)=='':
        messagebox.showwarning(title="Entrez", message=('No records found ...') ) 
    entrez_edit.insert('1.0',mess)
    
def Entrez0_Search():
    entrez_edit.delete('1.0', 'end')
    try:
        handle = Entrez.esearch(db=db0_Combo.get(), term=entrez_rec0_entry.get(), retmax=entrez_limit0_entry.get() )
    except BaseException as e:
        messagebox.showwarning(title="Entrez", message=str(e) )
        return False
    record = Entrez.read(handle)
    mess='Found:\n'
    for r in record["IdList"]:
        mess=mess+r+','
    entrez_edit.insert('1.0', mess)
    
def Swiss_Search():
    swiss_edit.delete('1.0', 'end')
    try:
        swiss_edit.delete('1.0', 'end')
        handle = ExPASy.get_sprot_raw(accessions_entry.get())
        swiss_edit.insert('1.0',handle.read())
    except BaseException as e:
        messagebox.showwarning(title="ExPASy", message=str(e) )
    
def Swiss0_Search():
    swiss_edit.delete('1.0', 'end')
    try:
        swiss_edit.delete('1.0', 'end')    
        handle = ExPASy.get_prosite_raw(prosite_entry.get())
        swiss_edit.insert('1.0',handle.read())
    except BaseException as e:
        messagebox.showwarning(title="Prosite", message=str(e) )
        
def Swiss1_Search():
    swiss_edit.delete('1.0', 'end')
    try:
        swiss_edit.delete('1.0', 'end')    
        handle = ScanProsite.scan(seq=prosite1_entry.get())
        swiss_edit.insert('1.0',handle.read())
    except BaseException as e:
        messagebox.showwarning(title="Prosite", message=str(e))
        
def TogoWS_Search():
    togows_edit.delete('1.0', 'end')
    label0_var.set(' Search in progress ...')
    window.update_idletasks()
    try:
        for id in TogoWS.search_iter(togo_Combo.get(), togo_rec_entry.get(), limit=int(togo_limit_entry.get())):
            handle=(TogoWS.entry(togo_Combo.get(), id, format=None, field=None))
            togows_edit.insert('1.0',handle.read())
    except BaseException as e:
        messagebox.showwarning(title="TogoWS", message=str(e) )
    label0_var.set('')
    
def Blast_Search():
    seq_file=seq_inuse_entry.get()
    xml_file=os.path.splitext(seq_inuse_entry.get())[0]+'.xml'     
    blast_edit.delete('1.0', 'end')
    label0_var.set(' Search in progress ...')
    window.update_idletasks()
    sequence_data = open(seq_file).read()
    #print(sequence_data)
    try:
        result_handle = NCBIWWW.qblast(blast_p_Combo.get(), blast_Combo.get(), sequence_data)
    except BaseException as e:
        messagebox.showwarning(title="Blast", message=str(e) )    
    with open(xml_file, 'w') as save_file: 
        blast_results = result_handle.read() 
        save_file.write(blast_results)
    E_VALUE_THRESH = float(blast_th_entry.get())
    mess=""
    for record in NCBIXML.parse(open(xml_file)): 
        if record.alignments:
            mess=mess+"Sequence: %s" % record.query[:130]+"\n"               
            for align in record.alignments: 
                for hsp in align.hsps: 
                    if hsp.expect < E_VALUE_THRESH: 
                        mess=mess+"   "+str(hsp.expect)+" %s " % align.title[:140]+"\n"                     
    blast_edit.insert('1.0',mess)  
    label0_var.set('')
    
def View_XML():
    xml_file=os.path.splitext(seq_inuse_entry.get())[0]+'.xml'
    #webbrowser.register('firefox',None,webbrowser.BackgroundBrowser(PATH_FIREFOX))
    #webbrowser.get('firefox').open('file:///'+xml_file)
    subprocess.run(["Notepad.exe", xml_file])
        
def Web(N): 
    if   (N=="0"):
      url = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    elif (N=="1"):
      url = 'http://www.ebi.ac.uk/Tools/msa/tcoffee/'  
    elif (N=="2"):
      url = 'http://www.ebi.ac.uk/Tools/msa/mafft/' 
    elif (N=="3"):
      url = 'http://www.ebi.ac.uk/goldman-srv/webprank/' 
    elif (N=="4"):
      url = 'http://toolkit.tuebingen.mpg.de/probcons' 
    elif (N=="5"):
      url = 'http://rna.informatik.uni-freiburg.de/LocARNA/Input.jsp'
    elif (N=="6"):
      url = 'http://rna.informatik.uni-freiburg.de/CARNA/Input.jsp'
    elif (N=="7"):
      url = 'http://rnacentral.org/'
      
    elif (N=="8"):
      url = 'http://www.phylogeny.fr/one_task.cgi?task_type=phyml/' 
    elif (N=="9"):
      url = 'http://www.phylogeny.fr/one_task.cgi?task_type=tnt' 
    elif (N=="10"):
      url = 'http://www.phylogeny.fr/one_task.cgi?task_type=bionj'
    elif (N=="11"):
      url = 'http://www.phylogeny.fr/one_task.cgi?task_type=mrbayes'
    elif (N=="12"):
      url = 'http://iubio.bio.indiana.edu/treeapp/treeprint-form.html'
      
    elif (N=="13"):
      url = 'http://primer3plus.com/cgi-bin/dev/primer3plus.cgi'
    elif (N=="14"):
      url = 'http://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi'
    elif (N=="15"):
      url = 'http://genome.ucsc.edu/cgi-bin/hgPcr'
    
    elif (N=="16"):
      url = 'http://mobyle.pasteur.fr/cgi-bin/portal.py#welcome'
    elif (N=="17"):
      url = 'https://usegalaxy.org/'
    elif (N=="18"):
      url = 'http://emboss.bioinformatics.nl/'
    elif (N=="19"):
      url = 'http://www.bioinformatics.org/sewer/'
    elif (N=="20"):
      url = 'http://www.bioservers.org/bioserver/'
    elif (N=="21"):
      url = 'http://www.trex.uqam.ca/index.php?action=home&tools=trex'
    elif (N=="22"):
      url = 'http://www.expasy.org/'
    elif (N=="23"):
      url = 'http://toolkit.lmb.uni-muenchen.de/'
    elif (N=="24"):
      url = 'http://www.lansys-sistemas.com/biolansys/sms2/'    
    webbrowser.register('firefox',None,webbrowser.BackgroundBrowser(PATH_FIREFOX))
    webbrowser.get('firefox').open(url)

def ExtApp(N): 
    if   (N=="0"):
        subprocess.run(["ExtApp\Clustalx.exe", ""])
    elif   (N=="1"):
        subprocess.run(["ExtApp\Aliview.exe", ""])
    elif   (N=="2"):
        subprocess.run(["ExtApp\FigTree.exe", ""])
        
def Seq_Ope(N):
    seq1_edit.delete('1.0', 'end')
    if   (N=="0"):      
        sequence = Seq(seq0_edit.get(1.0, tk.END))
        try:
         gc_seq=GC(sequence)
        except BaseException as e:
            messagebox.showwarning(title="GC content", message=str(e) )
        try:
            mw_seq=molecular_weight(sequence)
        except BaseException as e:
           messagebox.showwarning(title="Molecular weight", message=str(e) )
        mess="GC content: " + str(gc_seq) + " \n" + "Molecular Weight: " + str(mw_seq)
        seq1_edit.insert('1.0', mess)
        
    if   (N=="1"):       
        sequence = Seq(seq0_edit.get(1.0, tk.END))
        seq_list = list(sequence)
        seq_list.reverse()
        rev_seq = "".join(seq_list)
        seq1_edit.insert('1.0', rev_seq)        
    elif (N=="2"):       
        sequence = Seq(seq0_edit.get(1.0, tk.END))       
        seq_comp = sequence.complement()    
        seq1_edit.insert('1.0', seq_comp)
    elif   (N=="3"):       
        sequence = Seq(seq0_edit.get(1.0, tk.END))       
        seq_rna = transcribe(sequence)    
        seq1_edit.insert('1.0', seq_rna)
    elif   (N=="4"):
        try:         
           sequence = Seq(seq0_edit.get(1.0, tk.END).strip())      
           pro_seq = translate(sequence)  
           seq1_edit.insert('1.0', pro_seq)
        except BaseException as e:
           messagebox.showwarning(title="Translate Sequence", message=str(e) )

def Config(N):    
    if   (N=="0"):
        subprocess.run(["Notepad", "Config\clustalw.config"])
    elif (N=="1"):
        subprocess.run(["Notepad", "Config\clustalw.config"])
    elif (N=="2"):
        subprocess.run(["Notepad", "Config\clustalw.config"])
        
def Info(N):    
    if   (N=="0"):
        subprocess.run(["Notepad", "Docs\clustalw2.txt"])
    elif (N=="1"):
        subprocess.run(["Notepad", "Docs\clustalo.txt"])
    elif (N=="2"):
        subprocess.run(["Notepad", "Docs\muscle.txt"])
    elif (N=="3"):
        subprocess.run(["Notepad", "Docs\prank.txt"])
    elif (N=="4"):
        subprocess.run(["Notepad", "Docs\probcons.txt"])
    elif (N=="5"):
        subprocess.run(["Notepad", "Docs\dialign-tx.txt"])
    elif (N=="6"):
        subprocess.run(["Notepad", "Docs\pcma.txt"])
        
def Wiki():
    subprocess.run(["Notepad", "Docs\wiki.txt"])
        
# GUI
window = ThemedTk(theme='scidblue')
window.title("BioPython GUI V-0.3")
center_window(1240,600)
style = ttk.Style(window)
# Notebook style
s = ttk.Style()
s.configure('TNotebook.Tab', font=('Times','10' ) )
label_var = tk.StringVar()
label0_var = tk.StringVar()

# Main Frame
main_frame = ttk.Frame(window)
main_frame.pack(fill='both', expand=True)
# End Main Frame

# Top Frame
top_frame = tk.Frame(main_frame, bd=0, relief=tk.SOLID, padx=0, pady=0)
top_frame.pack(side=tk.TOP, fill=tk.X)

LoadSeq_btn = tk.Button(top_frame, width=8, height=1 ,text=' Load Seq ', font=('Lucida', 9,'bold'), command=Load_Seq)
CrearToolTip(LoadSeq_btn," Open Sequence File ")
LoadSeq_btn.grid(row=0, column=0, pady=3, padx=(2,2) )

SaveSeq_btn = tk.Button(top_frame, width=8, height=1 ,text=' Save Seq ', font=('Lucida', 9,'bold'), command=Save_Seq)
CrearToolTip(SaveSeq_btn," Save Sequence File ")
SaveSeq_btn.grid(row=0, column=1, pady=3, padx=(0,2) )

Align_Combo = ttk.Combobox(top_frame, state="readonly",width=8,font=('Lucida', 9,'bold') )
CrearToolTip(Align_Combo," Align Method ")
Align_Combo.set('Clustalw')
Align_Combo_Values=['Clustalw','ClustalO','Muscle','Prank','Probcons','Dialign','PCMA']
Align_Combo['values'] =  Align_Combo_Values 
Align_Combo.grid(sticky="W",row=0, column=2,   padx=(10,3), pady=(2,2))

AlignSeq_btn = tk.Button(top_frame, width=8, height=1 ,text=' Align ', font=('Lucida', 9,'bold'), command=Align_Seq)
CrearToolTip(AlignSeq_btn," Align Seq ")
AlignSeq_btn.grid(row=0, column=3, pady=3, padx=(0,5) )

Convert_Combo = ttk.Combobox(top_frame, state="readonly",width=12,font=('Lucida', 9,'bold') )
CrearToolTip(Convert_Combo," Convert Format ")
Convert_Combo.set('Fasta')
Convert_Combo_Values=['Fasta','Fasta2','Fastq','Fastq Sanger','Nexus','Genbank','embl','tab']
Convert_Combo['values'] =  Convert_Combo_Values 
Convert_Combo.grid(sticky="W",row=0, column=4,   padx=(10,3), pady=(2,2))

Convert_btn = tk.Button(top_frame, width=10, height=1 ,text=' Conversion ', font=('Lucida', 9,'bold'), command=Convert_Seq)
CrearToolTip(Convert_btn," Convert Seq ")
Convert_btn.grid(row=0, column=5, pady=3, padx=(0,5) )

process_label = tk.Label(top_frame, text="", textvariable=label0_var, bd=1,  anchor=tk.W, width=20,height=1,font=('Lucida', 9,'bold'),foreground = "red")
process_label.grid(row=0, column=6, pady=3, padx=(10,2) )
 
# End Top Frame

# Bottom Frame
bottom_frame = tk.Frame(main_frame,padx=0, pady=0)
bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

seq_inuse_label = tk.Label(bottom_frame,width=8,text="Seq File:",font=('Lucida', 9,'bold'),anchor=tk.W) 
seq_inuse_label.grid(row=0, column=0, pady=3, padx=(10,2) )
seq_inuse_entry = ttk.Entry(bottom_frame,width=75,font=('Lucida', 9,'bold'))   ; seq_inuse_entry.insert(0, "")     
seq_inuse_entry.grid(row=0, column=1, pady=3, padx=(2,2) )
format_label = tk.Label(bottom_frame, text="", textvariable=label_var, bd=1, relief=tk.SUNKEN, anchor=tk.W, width=10,height=1,font=('Lucida', 9,'bold'), background = "#7F8F9A",foreground = "white")
format_label.grid(row=0, column=2, pady=3, padx=(10,2) )

SeqType_Combo = ttk.Combobox(bottom_frame, state="readonly",width=8,font=('Lucida', 9,'bold') )
CrearToolTip(SeqType_Combo," Sequence Type ")
SeqType_Combo.set('DNA')
SeqType_Combo_Values=['DNA','PROTEIN']
SeqType_Combo['values'] =  SeqType_Combo_Values 
SeqType_Combo.grid(sticky="W",row=0, column=3,   padx=(10,3), pady=(2,2))
# End Bottom Frame

# Center Frame
Center_frame = tk.Frame(main_frame,width=180,bd=0, relief=tk.SOLID, padx=0, pady=0)
Center_frame.pack(side=tk.LEFT,fill="both", expand=1)
# Notebook
notebook0 = ttk.Notebook(Center_frame)
notebook0.pack(padx=0,pady=0,  expand = 1, fill ="both")
# create tabs
frame01 = ttk.Frame(notebook0) ; frame01.pack(fill='both', expand=True)
frame02 = ttk.Frame(notebook0) ; frame02.pack(fill='both', expand=True)
frame03 = ttk.Frame(notebook0) ; frame03.pack(fill='both', expand=True)
frame09 = ttk.Frame(notebook0) ; frame09.pack(fill='both', expand=True)
frame04 = ttk.Frame(notebook0) ; frame04.pack(fill='both', expand=True)
frame10 = ttk.Frame(notebook0) ; frame10.pack(fill='both', expand=True)
frame05 = ttk.Frame(notebook0) ; frame05.pack(fill='both', expand=True)
frame06 = ttk.Frame(notebook0) ; frame06.pack(fill='both', expand=True)
frame07 = ttk.Frame(notebook0) ; frame07.pack(fill='both', expand=True)
frame08 = ttk.Frame(notebook0) ; frame07.pack(fill='both', expand=True)
# add tabs to notebook
notebook0.add(frame01, text='  Sequences  ')
notebook0.add(frame02, text='  Align  ')
notebook0.add(frame03, text='  Phylo Tree  ')
notebook0.add(frame09, text='  Output  ')
notebook0.add(frame04, text='  Conversion  ')
notebook0.add(frame10, text='  Blastn  ')
notebook0.add(frame05, text='  Entrez  ')
notebook0.add(frame06, text='  SwissProt  ')
notebook0.add(frame07, text='  TogoWS  ')
 
notebook0.add(frame08, text='  Sequence  ')
# First tab
seq_edit = scrolledtext.ScrolledText(frame01,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
seq_edit.bind('<Button-3>',rClicker, add='')    
seq_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
seq_edit.message='seq_edit'
 
# Second tab
aln_edit = scrolledtext.ScrolledText(frame02,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
aln_edit.bind('<Button-3>',rClicker, add='')    
aln_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
aln_edit.message='aln_edit'

# PhyloTree tab
dnd_edit = scrolledtext.ScrolledText(frame03,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
dnd_edit.bind('<Button-3>',rClicker, add='')    
dnd_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
dnd_edit.message='dnd_edit'

# Output Tab
stdout_edit = scrolledtext.ScrolledText(frame09,font=('Lucida console', 10),  height=45,width=90,
           background = "#243642",foreground = "white",insertbackground="red")
stdout_edit.bind('<Button-3>',rClicker, add='')    
stdout_edit.pack(side=tk.LEFT,fill=tk.X,expand=True)
stdout_edit.message='stdout_edit'

stderr_edit = scrolledtext.ScrolledText(frame09,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
stderr_edit.bind('<Button-3>',rClicker, add='')    
stderr_edit.pack(side=tk.RIGHT,fill=tk.X,expand=True)
stderr_edit.message='stderr_edit'

# Fourth tab
conv_edit = scrolledtext.ScrolledText(frame04,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
conv_edit.bind('<Button-3>',rClicker, add='')    
conv_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
conv_edit.message='conv_edit'

# Five tab
entrez_top0_frame = tk.Frame(frame05,padx=0, pady=0)
entrez_top0_frame.pack(side=tk.TOP, fill=tk.X)

db0_Combo = ttk.Combobox(entrez_top0_frame, state="readonly",width=12,font=('Lucida', 9,'bold') )
CrearToolTip(db0_Combo," Database ")
db0_Combo.set('pubmed')
 
handle = Entrez.einfo()
record = Entrez.read(handle)
DBList=record["DbList"]
#print(DBList)
db0_Combo_Values=DBList 
db0_Combo['values'] =  db0_Combo_Values 
db0_Combo.grid(sticky="W",row=0, column=1,   padx=(10,3), pady=(2,2))

entrez_rec0_entry = ttk.Entry(entrez_top0_frame,width=60,font=('Lucida', 9,'bold'))   ; entrez_rec0_entry.insert(0, "diabetes+human")     
entrez_rec0_entry.grid(row=0, column=2, pady=3, padx=(2,2) )

entrez_limit0_entry = ttk.Entry(entrez_top0_frame,width=4,font=('Lucida', 9,'bold'))   ; entrez_limit0_entry.insert(0, "10")
entrez_limit0_entry.grid(row=0, column=3, pady=3, padx=(2,2) )
 
entrez0_btn = tk.Button(entrez_top0_frame, width=11, height=1 ,text=' Search Topic', font=('Lucida', 9,'bold'), command=Entrez0_Search)
CrearToolTip(entrez0_btn," Search Topic ")
entrez0_btn.grid(row=0, column=4, pady=3, padx=(2,2) )
###
entrez_top_frame = tk.Frame(frame05,padx=0, pady=0)
entrez_top_frame.pack(side=tk.TOP, fill=tk.X)

db_Combo = ttk.Combobox(entrez_top_frame, state="readonly",width=12,font=('Lucida', 9,'bold') )
CrearToolTip(db_Combo," Database ")
db_Combo.set('pubmed')
db_Combo['values'] =  db0_Combo_Values 
db_Combo.grid(sticky="W",row=0, column=1,   padx=(10,3), pady=(2,2))

id_entry = ttk.Entry(entrez_top_frame,width=80,font=('Lucida', 9,'bold'))   ; id_entry.insert(0, "6273291,6273290,6273289")     
id_entry.grid(row=0, column=2, pady=3, padx=(2,2) )
type_entry = ttk.Entry(entrez_top_frame,width=10,font=('Lucida', 9,'bold'))   ; type_entry.insert(0, "gb")     
type_entry.grid(row=0, column=3, pady=3, padx=(2,2) )
mode_entry = ttk.Entry(entrez_top_frame,width=10,font=('Lucida', 9,'bold'))   ; mode_entry.insert(0, "text")     
mode_entry.grid(row=0, column=4, pady=3, padx=(2,2) )
entrez_btn = tk.Button(entrez_top_frame, width=10, height=1 ,text=" Search ID's", font=('Lucida', 9,'bold'), command=Entrez_Search)
CrearToolTip(entrez_btn," Search ID's")
entrez_btn.grid(row=0, column=5, pady=3, padx=(2,2) )

entrez_edit = scrolledtext.ScrolledText(frame05,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
entrez_edit.bind('<Button-3>',rClicker, add='')    
entrez_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
entrez_edit.message='entrez_edit'

# Six tab
swiss_top_frame = tk.Frame(frame06,padx=0, pady=0)
swiss_top_frame.pack(side=tk.TOP, fill=tk.X)
swiss_top_frame0 = tk.Frame(frame06,padx=0, pady=0)
swiss_top_frame0.pack(side=tk.TOP, fill=tk.X)

swiss_label = tk.Label(swiss_top_frame,width=14,text="Accession/Name:",font=('Lucida', 9,'bold'),anchor=tk.W) 
swiss_label.grid(row=0, column=0, pady=3, padx=(10,2) )
accessions_entry = ttk.Entry(swiss_top_frame,width=15,font=('Lucida', 9,'bold'))   ; accessions_entry.insert(0, "O23729")     
accessions_entry.grid(row=0, column=1, pady=3, padx=(2,2) )
 
swiss_btn = tk.Button(swiss_top_frame, width=8, height=1 ,text=' Search ', font=('Lucida', 9,'bold'), command=Swiss_Search)
CrearToolTip(swiss_btn," Search SwissProt ")
swiss_btn.grid(row=0, column=5, pady=3, padx=(2,2) )

swiss0_label = tk.Label(swiss_top_frame,width=6,text=" Prosite:",font=('Lucida', 9,'bold'),anchor=tk.W) 
swiss0_label.grid(row=0, column=6, pady=3, padx=(10,2) )
prosite_entry = ttk.Entry(swiss_top_frame,width=15,font=('Lucida', 9,'bold'))   ; prosite_entry.insert(0, "PS50240")     
prosite_entry.grid(row=0, column=7, pady=3, padx=(2,2) )
swiss0_btn = tk.Button(swiss_top_frame, width=8, height=1 ,text=' Search ', font=('Lucida', 9,'bold'), command=Swiss0_Search)
CrearToolTip(swiss0_btn," Search Prosite ")
swiss0_btn.grid(row=0, column=8, pady=3, padx=(2,2) )

swiss1_label = tk.Label(swiss_top_frame0,width=8,text="Sequence:",font=('Lucida', 9,'bold'),anchor=tk.W) 
swiss1_label.grid(row=1, column=0, pady=3, padx=(10,2) )
prosite1_entry = ttk.Entry(swiss_top_frame0,width=105,font=('Lucida', 9,'bold'))   ;
prosite1_entry.insert(0, "MGSKRGISSRHHSLSSYEIMFAALFAILVVLCAGLIAVSCLTIKESQRGAALGQSHEARATFKITSGVTYNPNLQDKLSV")     
prosite1_entry.grid(row=1, column=1, pady=3, padx=(2,2) )
swiss1_btn = tk.Button(swiss_top_frame0, width=8, height=1 ,text=' Search ', font=('Lucida', 9,'bold'), command=Swiss1_Search)
CrearToolTip(swiss1_btn," Search Prosite ")
swiss1_btn.grid(row=1, column=2, pady=3, padx=(2,2) )

swiss_edit = scrolledtext.ScrolledText(frame06,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
swiss_edit.bind('<Button-3>',rClicker, add='')    
swiss_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
swiss_edit.message='swiss_edit'

# Togo WS tab
togows_top_frame = tk.Frame(frame07,padx=0, pady=0)
togows_top_frame.pack(side=tk.TOP, fill=tk.X)
togo_rec_label = tk.Label(togows_top_frame,width=7,text="TogoWS:",font=('Lucida', 9,'bold'),anchor=tk.W) 
togo_rec_label.grid(row=0, column=0, pady=3, padx=(10,2) )

togo_Combo = ttk.Combobox(togows_top_frame, state="readonly",width=12,font=('Lucida', 9,'bold') )
CrearToolTip(togo_Combo," Database ")
togo_Combo.set('pubmed')
togo_Combo_Values=['pubmed','nuccore','nucest','nucgss','nucleotide','protein','gene','omim','homologue','snp','mesh',
                   'embl','uniprot','uniparc','uniref100','uniref90','uniref50',
                   'ddbj','dad','pdb',
                   'compound','drug','enzyme','genes','glycan','orthology','reaction','module','pathway']
togo_Combo['values'] =  togo_Combo_Values 
togo_Combo.grid(sticky="W",row=0, column=1,   padx=(10,3), pady=(2,2))

togo_rec_entry = ttk.Entry(togows_top_frame,width=30,font=('Lucida', 9,'bold'))   ; togo_rec_entry.insert(0, "diabetes+human")     
togo_rec_entry.grid(row=0, column=2, pady=3, padx=(2,2) )

togo_limit_entry = ttk.Entry(togows_top_frame,width=4,font=('Lucida', 9,'bold'))   ; togo_limit_entry.insert(0, "10")
togo_limit_entry.grid(row=0, column=3, pady=3, padx=(2,2) )
 
togows_btn = tk.Button(togows_top_frame, width=8, height=1 ,text=' Search ', font=('Lucida', 9,'bold'), command=TogoWS_Search)
CrearToolTip(togows_btn," Open Sequence File ")
togows_btn.grid(row=0, column=4, pady=3, padx=(2,2) )

togows_edit = scrolledtext.ScrolledText(frame07,font=('Lucida console', 10),  height=45,
             background = "#243642",foreground = "white",insertbackground="red")
togows_edit.bind('<Button-3>',rClicker, add='')    
togows_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
togows_edit.message='togows_edit'

# Blast tab
blast_top_frame = tk.Frame(frame10,padx=0, pady=0)
blast_top_frame.pack(side=tk.TOP, fill=tk.X)
blast_label = tk.Label(blast_top_frame,width=10,text="Search Blast:",font=('Lucida', 9,'bold'),anchor=tk.W) 
blast_label.grid(row=0, column=0, pady=3, padx=(5,2) )

blast_p_Combo = ttk.Combobox(blast_top_frame, state="readonly",width=8,font=('Lucida', 9,'bold') )
CrearToolTip(blast_p_Combo," Program ")
blast_p_Combo.set('blastn')
blast_p_Combo_Values=['blastn','blastp', 'blastx', 'tblastn', 'tblastx']
blast_p_Combo['values'] =  blast_p_Combo_Values 
blast_p_Combo.grid(sticky="W",row=0, column=1,   padx=(10,3), pady=(2,2))

blast_Combo = ttk.Combobox(blast_top_frame, state="readonly",width=10,font=('Lucida', 9,'bold') )
CrearToolTip(blast_Combo," Database ")
blast_Combo.set('nt')
blast_Combo_Values=['nt','nr','RefSeq Select']
blast_Combo['values'] =  blast_Combo_Values 
blast_Combo.grid(sticky="W",row=0, column=2,   padx=(10,3), pady=(2,2))

blast_th_entry = ttk.Entry(blast_top_frame,width=8,font=('Lucida', 9,'bold'))   ; blast_th_entry.insert(0, "1e-20")     
blast_th_entry.grid(row=0, column=3, pady=3, padx=(2,2) )

blast_btn = tk.Button(blast_top_frame, width=8, height=1 ,text=' Search ', font=('Lucida', 9,'bold'), command=Blast_Search)
CrearToolTip(blast_btn," Blast  Search ")
blast_btn.grid(row=0, column=4, pady=3, padx=(2,2) )

blast0_btn = tk.Button(blast_top_frame, width=8, height=1 ,text=' View XML ', font=('Lucida', 9,'bold'), command=View_XML)
CrearToolTip(blast0_btn," View Blst XML result ")
blast0_btn.grid(row=0, column=5, pady=3, padx=(10,2) )

blast_edit = scrolledtext.ScrolledText(frame10,font=('Lucida console', 10),  height=45,
             background = "#243642",foreground = "white",insertbackground="red")
blast_edit.bind('<Button-3>',rClicker, add='')    
blast_edit.pack(side=tk.TOP,fill=tk.X,expand=True)
blast_edit.message='blast_edit'

# Sequence tab
seq0_top_frame = tk.Frame(frame08, bd=0, relief=tk.SOLID, padx=0, pady=0)
seq0_top_frame.pack(side=tk.TOP, fill=tk.X)
seq0_center_frame = tk.Frame(frame08, bd=0, relief=tk.SOLID, padx=0, pady=0)
seq0_center_frame.pack(side=tk.TOP, fill=tk.X)

bt00_seq0_btn = tk.Button(seq0_top_frame, width=8, height=1 ,text=' Stats ', font=('Lucida', 9,'bold'), command=lambda: Seq_Ope("0"))
CrearToolTip(bt00_seq0_btn," Sequence Stats ")
bt00_seq0_btn.grid(row=0, column=0, pady=3, padx=(2,2) )

bt01_seq0_btn = tk.Button(seq0_top_frame, width=8, height=1 ,text=' Reverse ', font=('Lucida', 9,'bold'), command=lambda: Seq_Ope("1"))
CrearToolTip(bt01_seq0_btn," Sequence Reverse ")
bt01_seq0_btn.grid(row=0, column=1, pady=3, padx=(2,2) )

bt02_seq0_btn = tk.Button(seq0_top_frame, width=10, height=1 ,text=' Complement ', font=('Lucida', 9,'bold'), command=lambda: Seq_Ope("2"))
CrearToolTip(bt02_seq0_btn," Sequence Complement ")
bt02_seq0_btn.grid(row=0, column=2, pady=3, padx=(2,2) )

bt03_seq0_btn = tk.Button(seq0_top_frame, width=12, height=1 ,text=' Transcription ', font=('Lucida', 9,'bold'), command=lambda: Seq_Ope("3"))
CrearToolTip(bt03_seq0_btn," DNA Sequence Transcription ")
bt03_seq0_btn.grid(row=0, column=3, pady=3, padx=(2,2) )

bt04_seq0_btn = tk.Button(seq0_top_frame, width=11, height=1 ,text=' Translation ', font=('Lucida', 9,'bold'), command=lambda: Seq_Ope("4"))
CrearToolTip(bt04_seq0_btn," DNA Sequence Translation ")
bt04_seq0_btn.grid(row=0, column=4, pady=3, padx=(2,2) )

seq0_edit = scrolledtext.ScrolledText(seq0_center_frame,font=('Lucida console', 10),  height=45, width=64,
           background = "#243642",foreground = "white",insertbackground="red")
seq0_edit.bind('<Button-3>',rClicker, add='')
seq0_edit.insert('1.0','AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAG')
seq0_edit.pack(side=tk.LEFT,fill=tk.X,expand=True)
seq0_edit.message='seq0_edit'

seq1_edit = scrolledtext.ScrolledText(seq0_center_frame,font=('Lucida console', 10),  height=45,
           background = "#243642",foreground = "white",insertbackground="red")
seq1_edit.bind('<Button-3>',rClicker, add='')    
seq1_edit.pack(side=tk.RIGHT,fill=tk.X,expand=True)
seq1_edit.message='seq1_edit'
 
# End Center Frame

# Menu bar
menubar = tk.Menu(window)

filemenu = tk.Menu(menubar, tearoff=0)
configmenu = tk.Menu(menubar, tearoff=0)
configmenu.add_command(label = "ClustalW", command=lambda: Config("0"))
configmenu.add_command(label = "ClustalO", command=lambda: Config("1"))
configmenu.add_command(label = "Muscle",   command=lambda: Config("2"))
configmenu.add_command(label = "Prank",    command=lambda: Config("3"))
configmenu.add_command(label = "Probcons", command=lambda: Config("4"))
filemenu.add_cascade(label="Config",menu=configmenu)
filemenu.add_separator()
filemenu.add_command(label="Exit"    , command=Close)
menubar.add_cascade(label="File"     , menu=filemenu)

alignmenu = tk.Menu(menubar, tearoff=0)
alignmenu.add_command(label="NCBIB: Blast",                 command=lambda: Web("0"))
alignmenu.add_command(label="EMBL-EBI: T-Coffee",          command=lambda: Web("1"))
alignmenu.add_command(label="EMBL-EBI: MAFFT",             command=lambda: Web("2"))
alignmenu.add_command(label="EMBL-EBI: webPRANK",          command=lambda: Web("3"))
alignmenu.add_command(label="MPlank: ProbCons ", command=lambda: Web("4"))
menubar.add_cascade(label="Web Align",     menu=alignmenu)

rnamenu = tk.Menu(menubar, tearoff=0)
rnamenu.add_command(label="LocARNA",     command=lambda: Web("5"))
rnamenu.add_command(label="CARNA",       command=lambda: Web("6"))
rnamenu.add_command(label="RNAcentral",  command=lambda: Web("7"))
menubar.add_cascade(label="Web RNA",     menu=rnamenu)

phylomenu = tk.Menu(menubar, tearoff=0)
phylomenu.add_command(label="PhyML",        command=lambda: Web("8"))
phylomenu.add_command(label="TNT",          command=lambda: Web("9"))
phylomenu.add_command(label="BioNJ",        command=lambda: Web("10"))
phylomenu.add_command(label="MrBayes",      command=lambda: Web("11"))
phylomenu.add_command(label="Phylodendron", command=lambda: Web("12"))
menubar.add_cascade(label="Web Phylo",      menu=phylomenu)

primermenu = tk.Menu(menubar, tearoff=0)
primermenu.add_command(label="PrimerPlus",       command=lambda: Web("13"))
primermenu.add_command(label="Primer-BLAST",         command=lambda: Web("14"))
primermenu.add_command(label="In-Silico PCR",       command=lambda: Web("15"))
menubar.add_cascade(label="Web Primer",     menu=primermenu)

suitesmenu = tk.Menu(menubar, tearoff=0)
suitesmenu.add_command(label="Mobyle",      command=lambda: Web("16"))
suitesmenu.add_command(label="Galaxy",      command=lambda: Web("17"))
suitesmenu.add_command(label="EMBOSS",      command=lambda: Web("18"))
suitesmenu.add_command(label="SeWeR",       command=lambda: Web("19"))
suitesmenu.add_command(label="BioServers",  command=lambda: Web("20"))
suitesmenu.add_command(label="T-REX",       command=lambda: Web("21"))
suitesmenu.add_command(label="ExPASy",      command=lambda: Web("22"))
suitesmenu.add_command(label="Bio/Toolkit", command=lambda: Web("23"))
suitesmenu.add_command(label="SMS2",        command=lambda: Web("24"))
menubar.add_cascade(label="Web Suites",     menu=suitesmenu)

extpmenu = tk.Menu(menubar, tearoff=0)
extpmenu.add_command(label="Clustalx",      command=lambda: ExtApp("0"))
extpmenu.add_command(label="AliView",       command=lambda: ExtApp("1"))
extpmenu.add_command(label="FigTree",       command=lambda: ExtApp("2"))
menubar.add_cascade(label="Ext Apps",       menu=extpmenu)

helpmenu = tk.Menu(menubar, tearoff=0)
infomenu = tk.Menu(menubar, tearoff=0)
infomenu.add_command(label = "ClustalW", command=lambda: Info("0"))
infomenu.add_command(label = "ClustalO", command=lambda: Info("1"))
infomenu.add_command(label = "Muscle",   command=lambda: Info("2"))
infomenu.add_command(label = "Prank",    command=lambda: Info("3"))
infomenu.add_command(label = "Probcons", command=lambda: Info("4"))
infomenu.add_command(label = "Dialign",  command=lambda: Info("5"))
infomenu.add_command(label = "PCMA",     command=lambda: Info("6"))
helpmenu.add_cascade(label="Info",menu=infomenu)
helpmenu.add_separator()
helpmenu.add_command(label="Wiki",     command=Wiki)
helpmenu.add_command(label="About...", command=About)
menubar.add_cascade(label="Help"     , menu=helpmenu)
window.config(menu=menubar)
# End Menu Bar

#window.state('zoomed') 
window.mainloop()

 