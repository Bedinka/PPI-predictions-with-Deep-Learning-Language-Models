let indexList = [];  // This array is assumed to be defined globally or in the same scope as checkIndex
//let exonNumDic = {};
let input_text = "";

window.indexList = indexList;

function dnaToAminoAcid(dnaSequence) {
    const codonToAA = {
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
    };

    let aminoAcids = "";
    for (let i = 0; i < dnaSequence.length; i += 3) {
        let codon = dnaSequence.substring(i, i + 3);
        let aminoAcid = codonToAA[codon] || '?'; // '?' as default for unknown codons
        aminoAcids += aminoAcid;
    }

    // Example of additional functionality as per your original Python script
    console.log(`DNA: ${dnaSequence}`);
    console.log(`Amino Acids: ${aminoAcids}`);
    
    // Additional functionality, including HTML generation, would go here

    return aminoAcids; // or HTML content, depending on your application's needs
}

function formatOutput(dnaSequence, aminoAcids, exonNumDic) {
    const lineLength = 20;  // amino acids per line
    let output = "";

    //console.log("$", dnaSequence);
    //console.log("$", aminoAcids);
    if ( indexList.length == 0){
        indexList.push(dnaSequence.length);
    }

    let foundDic = findIntronSiteFromAAs(aminoAcids); // Placeholder for your findIntronSiteFromAAs function

    for (let i = 0; i < aminoAcids.length; i += lineLength) {
        const aminoAcidsSegment = aminoAcids.substring(i, i + lineLength);
        const aminoAcidLine = aminoAcidsSegment.split('').join('  ');
        const dnaSegment = dnaSequence.substring(i * 3, i * 3 + lineLength * 3);
        //console.log("     " + aminoAcidLine);
        
        output += "      ";
        for (let j = 0; j < aminoAcidsSegment.length; j++) {
            //console.log(i + j);
            //console.log(aminoAcids.charAt(i + j));
            if (foundDic.includes(i + j)) {
                output += `<a href="/clicked" title="test"><font color="red">`;
                output += aminoAcids.charAt(i + j);
                output += "  ";
                output += `</font></a>`;
            } else {
                output += aminoAcids.charAt(i + j);
                output += "  ";
            }
        }
        
        output += "<br>";

        if (dnaSegment.length > 0) {
            //console.log("     " + dnaSegment);
            output += "      ";
            for (let j = 0; j < dnaSegment.length; j++) {
                const index = i * 3 + j;
                //console.log(index);
                //console.log(dnaSequence.charAt(index));

                // Check the exon number at the last index for evenness
                if (exonNumDic[index] % 2 === 0) {
                    output += "<font style='background-color: yellow'>";
                } else {
                    output += "<font style='background-color: pink'>";
                }
                if (dnaSequence.substring(index, index + 2) === "GG") {                
                    output += `<a href="javascript:htmx_click_site_test(${index});" title="Click to add/remove intron here - index:${index}" hx-post="/moclointron/htmx_click_site?index=${index}" hx-target="#exons" hx-swap="innerHTML"><font color="red">`;
                    output += dnaSequence.charAt(index);
                    output += `</font></a>`;
                } else {
                    //output += `<a href="javascript:htmx_click_site_test(${index});">`;
                    output += dnaSequence.charAt(index);
                    //output += `</a>`;
                }
                output += "</font>"
            }
            output += "<br>";
            // exon number
            output += "      ";
            for (let j = 0; j < dnaSegment.length; j++) {
                const index = i * 3 + j;
                output += exonNumDic[index] % 10;
            }
            output += "<br>";
            //console.log(`${i * 3 + 1} |` + '         |'.repeat(dnaSegment.length / 10)); //(dnaSegment.length / 30));
            output += `${String(i * 3 + 1).padStart(5)}|`;// + '         |' * 3;//(dnaSegment.length / 30);
            output += '         |'.repeat(dnaSegment.length / 10);
            output += "<br>";
        }
        //console.log("\n");
        output += "<br>";
    }
    return output;
}

function findIntronSiteFromAAs(cds) {
    let result = "";
    let foundAll = [];

    // JavaScript equivalent for Python's re.findall(), using matchAll() and a global regex.
    let found1 = Array.from(cds.matchAll(/.(?![C])(V|A|D|E|G)./g));
    let found2 = Array.from(cds.matchAll(/.(W|R|G)./g));
    let found3 = Array.from(cds.matchAll(/.G./g));

    // Extend foundAll with results from found1, found2, and found3
    foundAll = [...found1, ...found2, ...found3];

    let foundDic = {};
    foundAll.forEach(found => {
        let index = found.index;
        let value = found[0][1]; // Assuming the middle character is of interest
        //console.log(found);
        if (foundDic[index]) {
            foundDic[index].push(value);
        } else {
            foundDic[index] = [value];
        }
    });

    // Building the result string with marked elements
    for (let i = 0; i < cds.length; i++) {
        if (foundDic[i]) {
            result += `<a href="/clicked" title="test"><font color="red">`;
            result += cds[i];
            result += `</font></a>`;
        } else {
            result += cds[i];
        }
    }

    return result;  // Depending on your needs, you might return foundDic or result
}

// Example use
//let dnaSequence = "CCATCCCGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTAA";
//let aminoAcidSequence = dnaToAminoAcid(dnaSequence);
//console.log(`Translated Amino Acid Sequence: ${aminoAcidSequence}`);

//let output = formatOutput(dnaSequence, aminoAcidSequence);
//console.log(`Format output: ${output}`);

function setExonNum(indexList, mRNA) {
    let exonNum = 0;
    let exonNumDic = {};
    
    if ( indexList.length == 0){
        indexList.push(mRNA.length);
    }

    for (let i = 0; i < mRNA.length; i++) {
        exonNumDic[i] = exonNum;
        if (indexList.includes(i)) {
            exonNum += 1;
        }
    }

    return exonNumDic;
}

function update_Exon_Indexlist(newindexList){
    //alert(index);
    // Example use
    //let dnaSequence = "ATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTT";
    let dnaSequence = parent.mRNA.value;
    let aminoAcidSequence = dnaToAminoAcid(dnaSequence);
    console.log(`Translated Amino Acid Sequence: ${aminoAcidSequence}`);

    indexList = newindexList;
    console.log(`update_Exon_Indexlist: Index list: ${indexList}`);
    if ( indexList.length == 0){
        indexList.push(dnaSequence.length);
    }

    //checkIndex(index);
    update_Exon_Components(dnaSequence);
    exonNumDic = setExonNum(indexList, dnaSequence);

    let output = formatOutput(dnaSequence, aminoAcidSequence, exonNumDic);
    //output += `<br>from update_Exon_Indexlist: ${indexList}`;
    //console.log(`Format output: ${output}`);
    //alert(window.test.innerHTML);
    
    
    //window.test.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    //window.intron_output.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    //var scriptBlock = document.createElement('script');
    //scriptBlock.setAttribute("type","text/javascript");
    //scriptBlock.setAttribute("src", "/static/test.js");
    
    parent.frames['intron_output4'].document.body.innerHTML = "<pre style='font-size: 12pt;'>" + output + "</pre>";
    parent.window.indexList = indexList;
    parent.window.frames['intron_output4'].indexList = indexList;
    //window.frames['intron_output4'].document.head.appendChild(scriptBlock);
}

function htmx_click_site_test(index){
    //alert(index);
    // Example use
    //let dnaSequence = "ATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTT";
    let dnaSequence = parent.mRNA.value;
    let aminoAcidSequence = dnaToAminoAcid(dnaSequence);
    indexList = parent.window.indexList;
    console.log(`Translated Amino Acid Sequence: ${aminoAcidSequence}`);
    console.log(`htmx_click_site_test:Index list: ${indexList}`);

    if ( indexList.length == 0){
        indexList.push(dnaSequence.length);
    }

    checkIndex(index);
    update_Exon_Components(dnaSequence);
    exonNumDic = setExonNum(indexList, dnaSequence);

    let output = formatOutput(dnaSequence, aminoAcidSequence, exonNumDic);
    //output += `<br>from htmx_click_site_test: ${indexList}`;
    //console.log(`Format output: ${output}`);
    //alert(window.test.innerHTML);
    
    
    //window.test.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    //window.intron_output.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    //var scriptBlock = document.createElement('script');
    //scriptBlock.setAttribute("type","text/javascript");
    //scriptBlock.setAttribute("src", "/static/test.js");
    
    parent.frames['intron_output4'].document.body.innerHTML = "<pre style='font-size: 12pt;'>" + output + "</pre>";
    //console.log(`htmx_click_site_test:Index list: ${indexList}`);
    //window.frames['intron_output4'].document.head.appendChild(scriptBlock);
}

function update_Exon_Components(dnaSequence){
    let exon_element = parent.exons;
    //exon_element.innerHTML = "";
    let exon_text = "";
    let i = 0;
    let prev_index = 0;
    let first_exon_overhang_5 = parent.overhang_5.value;
    let last_exon_overhang_3 = parent.overhang_3.value;

    console.log(`update_Exon_Components: Index list: ${indexList}`);
    for (i = 0; i < indexList.length; i++) {
        //var input = document.createElement('input'); 
        //input.type = "text"; 
        //input.value = `Exon ${i}: ${indexList[i]}` 
        //exon_element.appendChild(input);
        //exon_text += `Exon ${i}:`;
        //exon_text += indexList[i];
        //exon_text += "<br>";
        let exon_len = indexList[i] - prev_index;
        let five_overhang = dnaSequence.substring(prev_index+1, prev_index+5);;
        let three_overhang = dnaSequence.substring(indexList[i]-3, indexList[i] + 1);;
        if (i % 2 === 0) {
            exon_text += "<font style='background-color: yellow'>";
        } else {
            exon_text += "<font style='background-color: pink'>";
        }
        exon_text += `&nbsp&nbsp&nbsp<input type="checkbox" class="checkbox_exons" id="exon_${i}" value="${indexList[i]}">&nbsp&nbsp <input type="text" value="Exon ${i}: ${prev_index+1} - ${indexList[i]} : (${exon_len})">`; // "<input type='text' value='test'><br>";
        exon_text += "</font>";
        //exon_text += `Length: ${exon_len}`; // "<input type='text' value='test'><br>";
        exon_text += "&nbsp&nbsp&nbsp&nbsp";
        if (i == 0){
            exon_text += `<input type="text" value="Overhang (5'): ${first_exon_overhang_5}">`; // "<input type='text' value='test'><br>";
        }
        else{
            exon_text += `<input type="text" value="Overhang (5'): ${five_overhang}">`; // "<input type='text' value='test'><br>";
        }
        if (i == indexList.length - 1){
            exon_text += `<input type="text" value="Overhang (3'): ${last_exon_overhang_3}">`; // "<input type='text' value='test'><br>";
        }
        else{
            exon_text += `<input type="text" value="Overhang (3'): ${three_overhang}">`; // "<input type='text' value='test'><br>";
        }
        exon_text += "<br>";
        prev_index = indexList[i];
    }
    //exon_text += `Exon ${i}:`;
    //exon_text += "${prev_index} - END";
    //exon_text += "<input type='checkbox'><br>";
    //exon_text += `<input type="checkbox"><input type="text" value="Exon ${i}: ${prev_index} - END">`; // "<input type='text' value='test'><br>";
    //exon_text += `Length: ${exon_len}<br>`; // "<input type='text' value='test'><br>";
    //exon_text += `<br>${indexList}`;
    exon_element.innerHTML = exon_text;
}

function checkIndex(index) {
    let foundIndex = indexList.indexOf(index);
    if (foundIndex !== -1) {
        indexList.splice(foundIndex, 1);  // Remove the index
    } else {
        indexList.push(index);  // Add the index
    }
    indexList.sort((a, b) => a - b);  // Sort the array numerically
    console.log(`Index list: ${indexList}`);
    return indexList;
}




//function submit_find_intron_site(index){
//    let dnaSequence = "ATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTT";
//    let aminoAcidSequence = dnaToAminoAcid(dnaSequence);
//    exonNumDic = setExonNum(dnaSequence);
//    let output = formatOutput(dnaSequence, aminoAcidSequence, exonNumDic);
//   window.intron_output.innerHTML = "<pre><h1>" + output + "</h1></pre>";
//}


function submit_find_intron_site(index) {
    input_text = window.mRNA.value;
    //alert(input_text);
    //let dnaSequence = "ATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTTATGTGTGCCGCTGCGGGCTCGGGATCCAGCGGGGGTGGCGGATTCTTCTTTGGAGGCGGCTTCTTT";
    let dnaSequence = input_text;
    let aminoAcidSequence = dnaToAminoAcid(dnaSequence);
    //console.log(`Translated Amino Acid Sequence: ${aminoAcidSequence}`);

    exonNumDic = setExonNum(indexList, dnaSequence);
    update_Exon_Components(dnaSequence);

    let output = formatOutput(dnaSequence, aminoAcidSequence, exonNumDic);
    //console.log(`Format output: ${output}`);
    //alert(window.test.innerHTML);
    //window.intron_output.innerHTML = "<iframe><pre><h1>" + output + "</h1></pre></iframe>";
    //window.intron_output4.document.body.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    var scriptBlock = document.createElement('script');
    scriptBlock.setAttribute("type","text/javascript");
    scriptBlock.setAttribute("src", "/static/test1.js");
    
    window.frames['intron_output4'].document.body.innerHTML = "<pre style='font-size: 12pt;'>" + output + "</pre>";
    window.frames['intron_output4'].document.head.appendChild(scriptBlock);
    //document.getElementsByTagName("intron_output4")[0].appendChild(scriptBlock);
    //window.intron_output.html("<pre><h1>" + output + "</h1></pre>");
    //window.intron_output3.innerHTML = "<pre><h1>" + output + "</h1></pre>";
    //window.frames['WDE_RTF'].document.body.innerHTML = "<pre>" + output + "</pre>";
}

function changeMoClo() {
    var x = document.getElementById("MoClo_select").value;
    //alert(x);
    if (x =="B3"){
        document.getElementById("overhang_5").value = "AATG";
        document.getElementById("overhang_3").value = "AGGT";
    }
    else if (x =="B3A"){
        document.getElementById("overhang_5").value = "AATG";
        document.getElementById("overhang_3").value = "ACGA";
    }
    else if (x =="B3B"){
        document.getElementById("overhang_5").value = "ACGA";
        document.getElementById("overhang_3").value = "AGGT";
    }
    else if (x =="B4"){
        document.getElementById("overhang_5").value = "AGGT";
        document.getElementById("overhang_3").value = "TTCG";
    }
    else if (x =="B5"){
        document.getElementById("overhang_5").value = "TTCG";
        document.getElementById("overhang_3").value = "GCTT";
    }
    else if (x =="B3-B4"){
        document.getElementById("overhang_5").value = "AATG";
        document.getElementById("overhang_3").value = "TTCG";
    }
    else if (x =="B3-B5"){
        document.getElementById("overhang_5").value = "AATG";
        document.getElementById("overhang_3").value = "GCTT";
    }
    else if (x =="B4-B5"){
        document.getElementById("overhang_5").value = "AGGT";
        document.getElementById("overhang_3").value = "GCTT";
    }
}

function reset_exons(){
    var newindexList = [];
    newindexList.push(indexList.pop());
    indexList = newindexList;
    window.indexList = newindexList;
    window.frames['intron_output4'].indexList = newindexList;
    console.log(`newindex list: ${newindexList}`);
    console.log(`Index list: ${window.frames['intron_output4'].indexList}`);
    update_Exon_Indexlist(newindexList);
}

function remove_selected_exons(){
    var newindexList = [];
    var checkedValue = null; 
    var inputElements = document.getElementsByClassName('checkbox_exons');
    console.log("remove_selected_exons");
    console.log(`Index list: ${indexList}`);
    console.log(`inputElements: ${inputElements.length}`);
    for(var i=0; i < inputElements.length-1; i++){
          if(inputElements[i].checked){
               //checkedValue = inputElements[i].value;
               console.log(`delete: ${i}, ${indexList[0]}, ${inputElements[i].value}`);
               //console.log(`newindex list: ${newindexList}`);
          }
          else{
               console.log(`no delete: ${i}, ${indexList[i]}, ${inputElements[i].value}`);
               newindexList.push(parseInt(inputElements[i].value));
               //console.log(`newindex list: ${newindexList}`);
          }
    }
    if(inputElements[i].checked){
        //checkedValue = inputElements[i].value;
        console.log(`delete: ${i}, ${indexList[0]}, ${inputElements[i].value}`);
        newindexList.pop();
        newindexList.push(parseInt(inputElements[i].value));
        //console.log(`newindex list: ${newindexList}`);
   }else{
        newindexList.push(parseInt(inputElements[i].value));
   }

    window.indexList = newindexList;
    window.frames['intron_output4'].indexList = newindexList;
    console.log(`newindex list: ${newindexList}`);
    console.log(`Index list: ${window.frames['intron_output4'].indexList}`);
    update_Exon_Indexlist(newindexList);
}

console.log("test1.js loaded 16!");



