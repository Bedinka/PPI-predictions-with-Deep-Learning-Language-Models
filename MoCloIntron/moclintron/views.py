from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader
from django.views.decorators.csrf import csrf_exempt
from django import forms

import re
from moclintron import codon_optimization_score

class Overhang_Form(forms.Form):
    overhangs = forms.CharField(widget=forms.Textarea)

class CDS_Form(forms.Form):
    cds = forms.CharField(widget=forms.Textarea)

class mRNA_Form(forms.Form):
    mRNA = forms.CharField(widget=forms.Textarea)

def getmRNA(cds):
    mRNA = codon_optimization_score.optimize_idt(cds)
    return mRNA

def search_all(pattern, text):
    result = []
    #for match in re.finditer(pattern, text):
    #    result.append(match.group())
    
    i = 0
    found = re.search(pattern, text)
    print(i, text, found)
    #for j in range(10): #while found != None:
    while True:
        if found == None:
            break
        i += (found.span()[0] + 1)
        result.append([i-1, found.group()])
        found = re.search(pattern, text[i:])
        print(i, text[i:], found)
        #result.append(found.group())

    return result

def wdeNumberToBase(n):
    return "TCAG"[n]

'''
function wdeNumberToBase(seq){
    var retSeq = "";
    switch (seq) {
        case 0: retSeq = "U";
            break;
        case 1: retSeq = "C";
            break;
        case 2: retSeq = "A";
            break;
        case 3: retSeq = "G";
            break;
    }
    return retSeq;
}'''

def wdeTranslateTripToAs(one, two, tre, transCode):

    return "M"
'''    var a = wdeBaseToNumber(one);
    var b = wdeBaseToNumber(two);
    var c = wdeBaseToNumber(tre);    
    var pos = a * 16 + b *4 + c;
    return wdeTranslate[bas][1].charAt(pos);'''

def wdeTranslateTripToStart(one, two, tre, transCode):
    return "M"

def wdeProteinOneThree(one):
    return "M"

def TranslationTable(mRNA, found_dic):
    ## setting
    wdeVTransCode = 1
    wdeVTransLetter = 1

    content = '<table border="0" style="line-height: 1.0; font-size: 80%;">'
    content += "<tr>\n<td></td><td></td><td colspan='4' style='text-align: center'>2nd Letter</td><td></td><td></td>\n</tr>\n"
    content += "<tr>\n<td></td><td></td><td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;U</td><td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C</td>"
    content += "<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A</td><td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G</td><td></td><td></td>\n</tr>\n"
    for k in range(4):
        for j in range(4):
            content += "<tr>\n"
            if ((k == 0) and (j == 0)):
	            content += "<td rowspan='16'>1st\nLetter&nbsp;&nbsp;&nbsp;</td>"

            if (((k == 0) or (k % 4)) and (j == 0)):
	            content += "<td rowspan='4'>&nbsp;" + wdeNumberToBase(k) + "&nbsp;&nbsp;&nbsp;</td>"
	        
            for i in range(4):
                one = wdeNumberToBase(k)
                two = wdeNumberToBase(i)
                tre = wdeNumberToBase(j)
                aas = wdeTranslateTripToAs(one, two, tre, wdeVTransCode)
                asStd = wdeTranslateTripToAs(one, two, tre, 0)
                asOrange = 0

                if (aas != asStd):
	                asOrange = 1
	            
                start = wdeTranslateTripToStart(one, two, tre, wdeVTransCode)

                if (start == "M"):
	                content += '<td>&nbsp;&nbsp;<span style="background-color:green">&nbsp;' + one + two + tre + "&nbsp;</span>"
                elif (aas == "*"):
	                content += '<td>&nbsp;&nbsp;<span style="background-color:red">&nbsp;' + one + two + tre + "&nbsp;</span>"
                else:
	                content += "<td>&nbsp;&nbsp;&nbsp;" + one + two + tre + "&nbsp;"
	            
                if (wdeVTransLetter == 1):
                    aas =  wdeProteinOneThree(aas)
	            
                if ((aas == "*") or (aas == "*  ")):
	                aas = "Stop"
	            
                if (asOrange):
	                content += '&nbsp;-&nbsp;<span style="background-color:orange">&nbsp;' + aas + "&nbsp;</span>&nbsp;&nbsp;</td>"
                else:
	                content += "&nbsp;-&nbsp;&nbsp;" + aas + "&nbsp;&nbsp;&nbsp;</td>"
	            
                if (i == 3):
                    content += "<td>&nbsp;&nbsp;&nbsp;" + wdeNumberToBase(j) + "&nbsp;</td>"
	            
            if ((k == 0) and (j == 0)):
	            content += "<td rowspan='16'>&nbsp;&nbsp;&nbsp;3rd\nLetter</td>"
	        
            content += "</tr>\n"

    content += "</table>"
    return content

def findIntronSiteFromAAs(cds):
    result = ""
    found_all = []

    found1 = search_all("\w[^C]V\w|\w[^C]A\w|\w[^C]D\w|\w[^C]E\w|\w[^C]G\w", cds)
    #result.append(["PATTERN 1","\w[^C]V\w|\w[^C]A\w|\w[^C]D\w|\w[^C]E\w|\w[^C]G\w",cds,found1])
    
    found2 = search_all("\wW\w|\wR\w|\wG\w", cds)
    #result.append(["PATTERN 2","\wW\w",cds,found2])

    found3 = search_all("\wG\w", cds)
    #result.append(["PATTERN 3","\wG\w",cds,found3])

    found_all.extend(found1)
    found_all.extend(found2)
    found_all.extend(found3)
    
    found_dic = {}
    for found in found_all:
        print(found)
        if found[0] in found_dic:
            found_dic[found[0]].append(found[1])
        found_dic[found[0]] = [found[1]]
        #result.append(found)
    
    for i in range(len(cds)):
        if i in found_dic:
            result += "<a href=\"/clicked\" title=\"test\"><font color=\"red\">"
            result += cds[i]
            result += "</font></a>"
        else:
            #result += "<font color=\"red\">"
            #result += mRNA[i]
            #result += "</font>"
            result += cds[i]

    return found_dic

def dna_to_amino_acid(dna_sequence):
    codon_to_aa = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    amino_acids = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        amino_acid = codon_to_aa.get(codon, '?')  # '?' as default for unknown codons
        amino_acids += amino_acid
    
    # Formatting output
    line_length = 20  # amino acids per line
    output = ""
    print("$", dna_sequence)
    print("$", amino_acids)
    found_dic = findIntronSiteFromAAs(amino_acids)

    for i in range(0, len(amino_acids), line_length):
        amino_acids_segment = amino_acids[i:i + line_length]
        amino_acid_line = "  ".join(amino_acids[i:i + line_length])
        dna_segment = dna_sequence[i * 3:(i + line_length) * 3]
        print("     " + amino_acid_line)
        
        output += "      "# + amino_acid_line
        for j in range(len(amino_acids_segment)):
            print(i+j)
            print(amino_acids[i+j])
            if i+j in found_dic:
                output += "<a href=\"/clicked\" title=\"test\"><font color=\"red\">"
                output += amino_acids[i+j]
                output += "  "
                output += "</font></a>"
            else:
                #result += "<font color=\"red\">"
                #result += mRNA[i]
                #result += "</font>"
                output += amino_acids[i+j]
                output += "  "
        
        output += "</br>"
        #print(f"  {i*3+1} |" + '         |' * 6 ) #(len(amino_acid_line)//10))
        #output += (f"  {i*3+1} |" + '         |' * 6 ) #(len(amino_acid_line)//10))
        #output += "</br>"

        if len(dna_segment) > 0:
            print("     " + dna_segment)

            output += "      "# + dna_segment
            for j in range(len(dna_segment)):
                print(i*3+j)
                index = i*3+j
                print(dna_sequence[i*3+j])
                if dna_sequence[i*3+j:i*3+j+2] == "GG":
                    # <a href="#" hx-post="/moclointron/htmx_click" onclick="alert('hi')" hx-target="#id_a" hx-swap="innerHTML">Click to swap content</a><br>
                    # output += "<a href=\"/clicked\" title=\"test\"><font color=\"red\">"
                    #output += ("<a href=\"javascript:htmx_click_site_test(%d);\" title=\"%d\" hx-post=\"/moclointron/htmx_click_site?index=%d\" hx-target=\"#exons\" hx-swap=\"innerHTML\"><font color=\"red\">" % (index,index,index))
                    output += ("<a href=\"javascript:htmx_click_site_test(%d);\"><font color=\"red\">" % index )
                    output += dna_sequence[i*3+j]
                    output += "</font></a>"
                else:
                    #output += ("<a href=\"javascript:htmx_click_site_test(%d);\">" % index )
                    output += dna_sequence[i*3+j]
                    #output += "</a>"
            output += "</br>"
            print(f" {i*3+1} |" + '         |' * (len(dna_segment)//10))
            output += "%5d" % (i*3+1)
            output += f"|" + '         |' * (len(dna_segment)//10)
            output += "</br>"
        print("\n")
        output += "</br>"
    return output
    
def htmx_click_site(request):
    index = request.GET.get('index')

    return HttpResponse("Click site %s" % index)

def htmx_mouse_entered(request):
    return HttpResponse("Mouse entered")

def htmx_messages(request):
    return HttpResponse("Message")

def htmx_click(request):
    return HttpResponse("Clicked")


def findIntronSite(mRNA):
    result = ""
    cds = mRNA#codon_optimization_score.translate(mRNA) #"1. ATG" + mRNA + "\n2. TAA"
    #found = re.findall("[^\wCV|^\wCA|!CD|!CE|!CG]", cds)
    #found = re.findall("[^C]A|[^C]G", cds)
    #found = re.search("[^C]A|[^C]G", cds)
    found_all = []

    found1 = search_all("\w[^C]V\w|\w[^C]A\w|\w[^C]D\w|\w[^C]E\w|\w[^C]G\w", cds)
    #result.append(["PATTERN 1","\w[^C]V\w|\w[^C]A\w|\w[^C]D\w|\w[^C]E\w|\w[^C]G\w",cds,found1])
    
    found2 = search_all("\wW\w|\wR\w|\wG\w", cds)
    #result.append(["PATTERN 2","\wW\w",cds,found2])

    found3 = search_all("\wG\w", cds)
    #result.append(["PATTERN 3","\wG\w",cds,found3])

    found_all.extend(found1)
    found_all.extend(found2)
    found_all.extend(found3)
    
    found_dic = {}
    for found in found_all:
        print(found)
        if found[0] in found_dic:
            found_dic[found[0]].append(found[1])
        found_dic[found[0]] = [found[1]]
        #result.append(found)
    
    for i in range(len(mRNA)):
        if i in found_dic:
            result += "<a href=\"/clicked\" title=\"test\"><font color=\"red\">"
            result += mRNA[i]
            result += "</font></a>"
        else:
            #result += "<font color=\"red\">"
            #result += mRNA[i]
            #result += "</font>"
            result += mRNA[i]
    return result

@csrf_exempt
def embed_test(request):
    template = loader.get_template("embed_test.html")
    print("embed_test")
    context = {
        "double_slide":...,
    }
    return HttpResponse(template.render(context, request))

@csrf_exempt
def description(request):
    template = loader.get_template("description.html")
    print("description")
    context = {
        "double_slide":...,
    }
    return HttpResponse(template.render(context, request))

@csrf_exempt
def double_slider(request):
    template = loader.get_template("double_slider.html")
    print("description")
    context = {
        "double_slide":...,
    }
    return HttpResponse(template.render(context, request))

# Create your views here.c:\Users\jaese\Dropbox\_CRAG\_Projects\2023_Moritz_Intron\codon_optimization_score.py
def example(request):
    return render(request, 'results.html')

def insert(request):
    print("request - insert")
    if request.method == 'POST':
        print("!====================")
        print(request.body)
        print("====================!")
        print(dir(request))
        print("request - insert - POST")
        if b"mRNA_output=" in request.body:
            form = CDS_Form(request.POST, request.FILES)
            if form.is_valid():
                print("True")
                print(form.cleaned_data)
                cds = form.cleaned_data['cds']
                #template = loader.get_template("index.html")
                context = {
                    "current_mRNA_output":getmRNA(cds),
                }
                return render(request, 'insert.html#mRNA_output', context) #HttpResponse(template.render(context, request))
        elif b"intron_output=" in request.body:
            form = mRNA_Form(request.POST, request.FILES)
            if form.is_valid():
                print("True")
                print(form.cleaned_data)
                mRNA = form.cleaned_data['mRNA']
                #template = loader.get_template("index.html")
                context = {
                    "current_intron_output":findIntronSite(mRNA),
                }
                html = render(request, 'insert.html#intron_output', context)
                print(dir(html))
                print(context["current_intron_output"])
                output = dna_to_amino_acid( mRNA)
                return HttpResponse('<div class="test" id="test"><pre><h1>%s</h1></pre></div>' % output )
                #return HttpResponse('<iframe name="intron_output2" id="intron_output2" style="border:#000000 1px solid; background-color: white; width:800px; height:500px;"><html><body><h1>Hello HttpResponse</h1></body></html></iframe>') #render(request, 'index.html#intron_output', context) #HttpResponse(template.render(context, request))
                return HttpResponse('<div class="test"><h1>%s</h1></div>' % context["current_intron_output"]) #render(request, 'index.html#intron_output', context) #HttpResponse(template.render(context, request))
        else:
            form = Overhang_Form(request.POST, request.FILES)
            if form.is_valid():
                current_overhangs = form.cleaned_data['overhangs']
                fidelity_output = test(current_overhangs)
                print(fidelity_output)
                context = {
                    "current_overhangs": current_overhangs, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
                    "current_fidelity_output": fidelity_output, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
                }
                return render(request, 'insert.html#fidelity_output', context) #HttpResponse(template.render(context, request))

    template = loader.get_template("insert.html")
    context = {
        "latest_question_list": ...,
    }
    return HttpResponse(template.render(context, request))


# Create your views here.c:\Users\jaese\Dropbox\_CRAG\_Projects\2023_Moritz_Intron\codon_optimization_score.py
def index(request):
    print("request - index")
    if request.method == 'POST':
        print("!====================")
        print(request.body)
        print("====================!")
        print(dir(request))
        print("request - index - POST")
        if b"mRNA_output=" in request.body:
            form = CDS_Form(request.POST, request.FILES)
            if form.is_valid():
                print("True")
                print(form.cleaned_data)
                cds = form.cleaned_data['cds']
                #template = loader.get_template("index.html")
                context = {
                    "current_mRNA_output":getmRNA(cds),
                }
                return render(request, 'index.html#mRNA_output', context) #HttpResponse(template.render(context, request))
        elif b"intron_output=" in request.body:
            form = mRNA_Form(request.POST, request.FILES)
            if form.is_valid():
                print("True")
                print(form.cleaned_data)
                mRNA = form.cleaned_data['mRNA']
                #template = loader.get_template("index.html")
                context = {
                    "current_intron_output":findIntronSite(mRNA),
                }
                html = render(request, 'index.html#intron_output', context)
                print(dir(html))
                print(context["current_intron_output"])
                output = dna_to_amino_acid( mRNA)
                return HttpResponse('<div class="test" id="test"><pre><h1>%s</h1></pre></div>' % output )
                #return HttpResponse('<iframe name="intron_output2" id="intron_output2" style="border:#000000 1px solid; background-color: white; width:800px; height:500px;"><html><body><h1>Hello HttpResponse</h1></body></html></iframe>') #render(request, 'index.html#intron_output', context) #HttpResponse(template.render(context, request))
                return HttpResponse('<div class="test"><h1>%s</h1></div>' % context["current_intron_output"]) #render(request, 'index.html#intron_output', context) #HttpResponse(template.render(context, request))
        else:
            form = Overhang_Form(request.POST, request.FILES)
            if form.is_valid():
                current_overhangs = form.cleaned_data['overhangs']
                fidelity_output = test(current_overhangs)
                print(fidelity_output)
                context = {
                    "current_overhangs": current_overhangs, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
                    "current_fidelity_output": fidelity_output, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
                }
                return render(request, 'index.html#fidelity_output', context) #HttpResponse(template.render(context, request))

    template = loader.get_template("index.html")
    context = {
        "latest_question_list": ...,
    }
    return HttpResponse(template.render(context, request))

@csrf_exempt
def results(request):
    template = loader.get_template("results.html")
    form = Overhang_Form(request.POST, request.FILES)
    if form.is_valid():
        print("True")
        print(form.cleaned_data)
        current_overhangs = form.cleaned_data['overhangs']
        output = test(current_overhangs)
        context = {
            "current_overhangs": current_overhangs, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
            "current_output": output, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
        }
    else:
        print("False")
        context = {
            "current_overhangs": "",
            "current_output": "",
        }
    return HttpResponse(template.render(context, request))


def optimize_codons(cds):
    mRNA = ""

    return mRNA

@csrf_exempt
def optimize_codons(request):
    template = loader.get_template("results.html")
    form = Overhang_Form(request.POST, request.FILES)
    if form.is_valid():
        print("optimize_codons")
        print(form.cleaned_data)
        current_overhangs = form.cleaned_data['overhangs']
        output = test(current_overhangs)
        context = {
            "current_overhangs": current_overhangs, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
            "current_output": output, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
        }
    else:
        print("False")
        context = {
            "current_overhangs": "",
            "current_output": "",
        }
    return HttpResponse(template.render(context, request))

import json


@csrf_exempt
def calculate_fidelity_axios(request):
    template = loader.get_template("results.html")
    print(request.body)
    print(request.POST.get('data'))
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    overhangs = body['overhangs']
    print(overhangs)
    fidelity_output = calc_fidelity_from_overhangs(overhangs)
    print(fidelity_output)
    return HttpResponse(fidelity_output)

@csrf_exempt
def optimize_codons_axios(request):
    template = loader.get_template("results.html")
    print(request.body)
    print(request.POST.get('data'))
    print(dir(request))
    body_unicode = request.body.decode('utf-8')
    body = json.loads(body_unicode)
    protein = body['protein']
    mRNA = getmRNA(protein)
    context={
        "current_overhangs": "current_overhangs",
        "current_output": "current_output"
    }
    '''
    form = Overhang_Form(request.POST, request.FILES)
    if form.is_valid():
        print("optimize_codons")
        print(form.cleaned_data)
        current_overhangs = form.cleaned_data['overhangs']
        output = test(current_overhangs)
        context = {
            "current_overhangs": current_overhangs, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
            "current_output": output, #form.cleaned_data['overhangs'], #request.POST["overhangs"],
        }
    else:
        print("False")
        context = {
            "current_overhangs": "",
            "current_output": "",
        }
    '''
    #mRNA = "ATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTT"
    return HttpResponse(mRNA)
    return HttpResponse(template.render(context, request))

def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

def specificity_with( tableX, overhang, overhang_list ):
    cnt = 0
    total = 0
    rc = rev_compl(overhang)
    cnt = tableX[(overhang, rc)]
    total = tableX[(overhang, rc)]
    for overhang2 in overhang_list:
        if overhang2 == rc: continue
        total += tableX[(overhang, overhang2)]
    return cnt, total, cnt/float(total)

def get_specificity_list(table, overhang_set):
    all_overhangs = set()
    rc_overhangs = set()
    for overhang1 in overhang_set:
        #ic(overhang1)
        #ic(IUPACT2Seqs(overhang1))
        rc = rev_compl(overhang1)
        all_overhangs.add(overhang1)
        rc_overhangs.add(rc)
    
    #ic(all_overhangs)
    #ic(rc_overhangs)
    
    all_rc_overhangs = set()
    all_rc_overhangs.update( all_overhangs )
    all_rc_overhangs.update( rc_overhangs )

    #ic(all_rc_overhangs)

    specificity_list = []
    
    for overhang in all_overhangs:
        [cnt, total, freq ] = specificity_with(table, overhang, all_rc_overhangs)
        specificity_list.append( [freq, overhang] )
        #ic(overhang, cnt, total, freq)
    
    specificity_list.sort(reverse=True)
    return specificity_list

def specificity_within_group( tableX, Input_overhangs, overhangs, return_fidelity_only = False ):    
    output = ""                 
    unspecific_total = 0
    specific_total = 0
    fidelity = 1.0
    for overhang in overhangs:
        specific = 0
        unspecific = 0
        for overhang2 in overhangs:
            #overhang2 = rev_compl(overhang2)
            if overhang == rev_compl(overhang2): #overhang2: #rev_compl(overhang2):
                if overhang == rev_compl(overhang):
                    print("** Specific ** ", overhang, overhang2, tableX[(overhang, overhang2)] )
                    output += "** Specific ** " + overhang + " " + overhang2 + " " + str(tableX[(overhang, overhang2)]) + "\n"
                    specific += tableX[(overhang, overhang2)]
                    unspecific += tableX[(overhang, overhang2)]
                else:
                    print("** Specific ** ", overhang, overhang2, tableX[(overhang, overhang2)] )
                    output += "** Specific ** " + overhang + " " + overhang2 + " " + str(tableX[(overhang, overhang2)]) + "\n"
                    specific += tableX[(overhang, overhang2)]
            else:
                print("Unspecific ", overhang, overhang2, tableX[(overhang, overhang2)] )
                output += "Unspecific " + overhang + " " + overhang2 + " " + str(tableX[(overhang, overhang2)]) + "\n"
                unspecific += tableX[(overhang, overhang2)]
        print("## Specificity =", specific, unspecific, specific/float(specific+unspecific))
        output += "## Specificity = " + str(specific) + " " + str(unspecific) + " " + str(specific/float(specific+unspecific)) + "\n"
        if overhang in Input_overhangs:
            fidelity *= specific/float(specific+unspecific)
        print("\n")
        output += "\n"
        
        unspecific_total += unspecific
        specific_total += specific

    total = specific_total + unspecific_total
    print("## Total specificity =", overhangs, specific_total, unspecific_total)
    print("## Total specificity =", specific_total, unspecific_total, specific_total/float(total))
    output += "## Total specificity = " + str(specific_total) + " " + str(unspecific_total) + " " + str(specific_total/float(total)) + "\n"
    output += "## Fidelity = " + str(fidelity) + "\n"
    
    '''
    rc = rev_compl(overhang)
    nt_set = set( all_nt(len(overhang)) )
    nt_set.discard(rc)

    cnt = 0
    for overhang2 in nt_set:
        cnt += tableX[(overhang, overhang2)]
    
    total = tableX[(overhang, rc)] + cnt 
    '''
    
    if return_fidelity_only:
        return fidelity
    
    return output #specific_total, unspecific_total, specific_total/float(total)

def read_ligation_freq_file(filepath, token="\t"):
    table = {}
    f = open(filepath)
    columns = f.readline()[:-1].split(token)[1:]
    for line in f:
        fields = line[:-1].split(token)
        for i in range(len(columns)):
            #table[(fields[0],rev_compl(columns[i]))] = float(fields[i+1])
            table[(fields[0], columns[i])] = float(fields[i+1])
    return table

#BsaI_HFv2 = read_ligation_freq_file("C:\\Users\\jaese\\Dropbox\\Code3\\MoCloIntron\\myproject\\moclintron\\BsaI-HFv2.txt") # 4NTs overhang
#print("read BsaI-HFv2.txt")
#GoldenGateLigationInfo = read_ligation_freq_file("C:\\Users\\jaese\\Dropbox\\Code3\\MoCloIntron\\myproject\\moclintron\\BsaI-HFv2.txt") # 4NTs overhang
#print("read BsaI-HFv2.txt")

#GoldenGateLigationInfo = read_ligation_freq_file("C:\\Users\\jaese\\Dropbox\\Code3\\MoCloIntron\\myproject\\moclintron\\BbsI-HF.txt") # 4NTs overhang
GoldenGateLigationInfo = read_ligation_freq_file("./moclintron/BbsI-HF.txt") # 4NTs overhang
print("read BbsI-HF.txt")
print(GoldenGateLigationInfo[("ACGG", "CCGC")])
print(GoldenGateLigationInfo[("CCGC", "ACGG")])
print(GoldenGateLigationInfo[("CCGG", "TCGG")])
print(GoldenGateLigationInfo[("TCGG", "CCGG")])
print(GoldenGateLigationInfo[("TACC", "TACC")])
print(GoldenGateLigationInfo[("TACC", "GGTT")])
print(GoldenGateLigationInfo[("TACC", "GGTA")])

#GoldenGateLigationInfo2 = read_ligation_freq_file("C:\\Users\\jaese\\Dropbox\\Code3\\MoCloIntron\\myproject\\moclintron\\BsaI-HFv2.txt") # 4NTs overhang
GoldenGateLigationInfo2 = read_ligation_freq_file("./moclintron/BsaI-HFv2.txt") # 4NTs overhang
print("read BsaI-HFv2.txt")
print(GoldenGateLigationInfo2[("ACGG", "CCGC")])
print(GoldenGateLigationInfo2[("CCGC", "ACGG")])
print(GoldenGateLigationInfo2[("CCGG", "TCGG")])
print(GoldenGateLigationInfo2[("TCGG", "CCGG")])
print(GoldenGateLigationInfo2[("TACC", "TACC")])
print(GoldenGateLigationInfo2[("TACC", "GGTT")])
print(GoldenGateLigationInfo2[("TACC", "GGTA")])

GoldenGateLigationInfo3 = read_ligation_freq_file("./moclintron/sb8b00333_si_002/Supplemental Data/FileS08_T7_18h_37C.csv", token=",") # 4NTs overhang
print("read FileS08_T7_18h_37C.csv")
print(GoldenGateLigationInfo3[("ACGG", "CCGC")])
print(GoldenGateLigationInfo3[("CCGC", "ACGG")])
print(GoldenGateLigationInfo3[("CCGG", "TCGG")])
print(GoldenGateLigationInfo3[("TCGG", "CCGG")])
print(GoldenGateLigationInfo3[("TACC", "TACC")])
print(GoldenGateLigationInfo3[("TACC", "GGTT")])
print(GoldenGateLigationInfo3[("TACC", "GGTA")])

GoldenGateLigationInfo4 = read_ligation_freq_file("./moclintron/sb8b00333_si_002/Supplemental Data/FileS06_T7_18h_25C.csv", token=",") # 4NTs overhang
print("read FileS06_T7_18h_25C.csv")
print(GoldenGateLigationInfo4[("ACGG", "CCGC")])
print(GoldenGateLigationInfo4[("CCGC", "ACGG")])
print(GoldenGateLigationInfo4[("CCGG", "TCGG")])
print(GoldenGateLigationInfo4[("TCGG", "CCGG")])
print(GoldenGateLigationInfo4[("TACC", "TACC")])
print(GoldenGateLigationInfo4[("TACC", "GGTT")])
print(GoldenGateLigationInfo4[("TACC", "GGTA")])




def test(overhangs):
    output = "==============================\n"
    output += "1. Add Overhang Sequences\n"
    output += "==============================\n\n"

    Input_overhangs = set()
    OvoA_intron_overhang_set = set()
    overhangs = overhangs.upper()

    for overhang1 in overhangs.split("\n"):
        if "," in overhang1:
            for overhang in overhang1.split(","):
                overhang = overhang.strip()
                if len(overhang) != 4: continue
                Input_overhangs.add(overhang)
        else:
            overhang = overhang1.strip()
            if len(overhang) != 4: continue
            Input_overhangs.add(overhang)
        
    for overhang in Input_overhangs:
        print(overhang)
        output += "# Add Overhang Sequence:" + overhang + "\n"
        OvoA_intron_overhang_set.add(overhang)
        output += "# Add Reverse Complement Sequence of Overhang: " + rev_compl(overhang) + "\n"
        OvoA_intron_overhang_set.add(rev_compl(overhang))
        output += "\n"
        
        #else:
        #    print(overhang1)
        #    overhang = overhang1.strip()
        #    if len(overhang) != 4: continue
        #    Input_overhangs.add(overhang)
        #    output += "# Add Overhang Sequence:" + overhang + "\n"
        #    OvoA_intron_overhang_set.add(overhang)
        #    output += "# Add Reverse Complement Sequence of Overhang: " + rev_compl(overhang) + "\n"
        #    OvoA_intron_overhang_set.add(rev_compl(overhang))
        #    output += "\n"
    '''
    OvoA_intron_overhang_set = set()
    OvoA_intron_overhang_set.add("GACG")
    #OvoA_intron_overhang_set.add(rev_compl("GGCG"))    # ## Specificity = 280.0 97.0 0.7427055702917772
    OvoA_intron_overhang_set.add(rev_compl("GGAG"))
    
    OvoA_intron_overhang_set.add("CAAG")
    OvoA_intron_overhang_set.add(rev_compl("GCCA"))
    
    OvoA_intron_overhang_set.add("CGCG")
    OvoA_intron_overhang_set.add(rev_compl("GCGG"))
    
    OvoA_intron_overhang_set.add("TACG")
    OvoA_intron_overhang_set.add(rev_compl("GCCC"))
    
    OvoA_intron_overhang_set.add("CGAG")
    OvoA_intron_overhang_set.add(rev_compl("GTGG"))
    '''
    Final_overhangs = set(OvoA_intron_overhang_set)

    output += "==============================\n"
    output += "2. Simulating ligation events\n"
    output += "==============================\n\n"
    
    #for overhang in OvoA_intron_overhang_set:
    #    Final_overhangs.add(rev_compl(overhang))

    #specificity_list = get_specificity_list( GoldenGateLigationInfo, Final_overhangs)
    #print(specificity_list)
    #print( specificity_within_group(GoldenGateLigationInfo, Final_overhangs) )

    output += specificity_within_group(GoldenGateLigationInfo, Input_overhangs, Final_overhangs)
    return output



def calc_fidelity_from_overhangs(overhangs):
    output = "==============================\n"
    output += "1. Add Overhang Sequences\n"
    output += "==============================\n\n"

    Input_overhangs = set()
    OvoA_intron_overhang_set = set()
    
    for overhang1 in overhangs:
        if "," in overhang1:
            for overhang in overhang1.split(","):
                overhang = overhang.strip()
                if len(overhang) != 4: continue
                Input_overhangs.add(overhang)
        else:
            overhang = overhang1.strip()
            if len(overhang) != 4: continue
            Input_overhangs.add(overhang)
        
    for overhang in Input_overhangs:
        print(overhang)
        output += "# Add Overhang Sequence:" + overhang + "\n"
        OvoA_intron_overhang_set.add(overhang)
        output += "# Add Reverse Complement Sequence of Overhang: " + rev_compl(overhang) + "\n"
        OvoA_intron_overhang_set.add(rev_compl(overhang))
        output += "\n"
        
        #else:
        #    print(overhang1)
        #    overhang = overhang1.strip()
        #    if len(overhang) != 4: continue
        #    Input_overhangs.add(overhang)
        #    output += "# Add Overhang Sequence:" + overhang + "\n"
        #    OvoA_intron_overhang_set.add(overhang)
        #    output += "# Add Reverse Complement Sequence of Overhang: " + rev_compl(overhang) + "\n"
        #    OvoA_intron_overhang_set.add(rev_compl(overhang))
        #    output += "\n"

    Final_overhangs = set(OvoA_intron_overhang_set)

    output += "==============================\n"
    output += "2. Simulating ligation events\n"
    output += "==============================\n\n"
    
    print(output)
    #for overhang in OvoA_intron_overhang_set:
    #    Final_overhangs.add(rev_compl(overhang))

    #specificity_list = get_specificity_list( GoldenGateLigationInfo, Final_overhangs)
    #print(specificity_list)
    #print( specificity_within_group(GoldenGateLigationInfo, Final_overhangs) )

    fidelity = specificity_within_group(GoldenGateLigationInfo, Input_overhangs, Final_overhangs, True)
    print("Fidelity: {:0.3f}".format(fidelity))
    return "{:0.3f}".format(fidelity)