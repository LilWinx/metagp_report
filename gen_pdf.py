from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib import utils, colors
from reportlab.platypus import Table, TableStyle

def pdf_template(outdir, in_dict, used_reference):
    outpath = outdir + "/report.pdf"
    c = canvas.Canvas(outpath, pagesize=A4)

    # Letter head stuff
    logo_img = utils.ImageReader(in_dict.get("py_logo_ph", ""))
    metagp_img = utils.ImageReader(in_dict.get("py_metagplogo_ph", ""))

    c.drawImage(logo_img, 10, A4[1] - 54, width=110, height=44, mask="auto")
    c.drawImage(metagp_img, (A4[0] / 2) - 55, A4[1] - 54, width=110, height=30, mask="auto")

    address = [
        "Centre for Infectious Diseases and Microbiology Laboratory Services",
        "Level 3, NSW Health Pathology, ICPMR, Westmead Hospital",
        "Westmead, NSW, 2145, Australia",
    ]
    address_y = A4[1] - 30
    c.setFont("Helvetica", 5)
    for line in address:
        c.drawRightString(A4[0] - 10, address_y, line)
        address_y -= 6

    # Title
    title = "Metagenomics Report"
    c.setFont("Helvetica", 12)
    c.drawCentredString(A4[0] / 2, A4[1] - 80, title)
    c.translate(0, -60) # add space after
    
    # Patient Info table
    pt_metadata = [['Patient Name:', 'py_pn_ph', 'MRN:', 'py_mrn_ph'],
                    ["Accession:", "py_acc_ph", "Date of Collection:", "py_doc_ph"],
                    ["Date of Birth:", "py_dob_ph", "Request Date:", "py_reqdate_ph"],
                    ["Sequence Date:","py_seqdate_ph","Run Name:","py_runname_ph"],
                    ["Date of Report:","py_repdate_ph", "WGS ID:", "py_wgsid_ph"]]
    for row_idx, row in enumerate(pt_metadata):
        for col_idx, cell in enumerate(row):
            if cell in in_dict:
                pt_metadata[row_idx][col_idx] = in_dict[cell]

    c.setFont("Helvetica", 9)
    table = Table(pt_metadata, rowHeights=4*3.5)
    
    style = TableStyle([
        ('FONTNAME', (0, 0), (4, 4), 'Helvetica'),
        ('FONTSIZE', (0, 0), (4, 4), 8),
        ('RIGHTPADDING', (1, 0), (1, -1), 40),
        ('BOTTOMPADDING', (0, 0), (4, 4), 1),
    ])
    
    table.setStyle(style)
    table.wrapOn(c, 0,0)
    table.drawOn(c, A4[0] / 2 - 175, A4[1] - 110)
    
    # draw dividing line
    x_line, y_line = 50, A4[1] - 130
    x_line2, y_line2 = 550, A4[1] - 130
    c.line(x_line, y_line, x_line2, y_line2)

    # Results section  
    clinical_results_list = [
        "Date of Illness Onset: 'py_date_onset_ph'",
        "Suspected Pathogen: 'py_suspathogen_ph'",
        "Sample Type: 'py_sampletype_ph'",
        "Clinician Comments: 'py_clincomm_ph'",
    ]

    meta_summary_list = [
        "Interpretation: 'py_riinterpret_ph'",
        "Reason: 'py_rireason_ph'"
    ]
    line_count = 0
    header = "Clinical Request Information:"
    header2 = "Metagenomics Summary:"
    header3 = "Methods:"
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, A4[1]- 150, header)
    text_width = c.stringWidth(header, "Helvetica-Bold", 8)
    c.line(50, A4[1] - 152, 50 + text_width, A4[1] - 152)
    
    c.setFont("Helvetica", 8)

    # Your clinical results
    clinical_buffer = A4[1] - 170
    for result in clinical_results_list:
        placeholder = result.split(": ")[1].strip("'")
        if placeholder in in_dict:
            result = result.replace(f"'{placeholder}'", in_dict[placeholder])
            numlines = draw_wrapped_text(c, result, clinical_buffer, 8)
            clinical_buffer -= 15
            line_count += numlines
    
    # meta summary
    summary_buffer = clinical_buffer - (line_count * 8)
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, summary_buffer, header2)
    text_width2 = c.stringWidth(header2, "Helvetica-Bold", 8)
    c.line(50, summary_buffer - 2, 50 + text_width2, summary_buffer - 2)
    line_count = 0
    
    summary_buffer = summary_buffer - 15
    for summary in meta_summary_list:
        placeholder = summary.split(": ")[1].strip("'")
        if placeholder in in_dict:
            summary = summary.replace(f"'{placeholder}'", in_dict[placeholder])
            numlines = draw_wrapped_text(c, summary, summary_buffer, 8)
            summary_buffer -= 15
            line_count += numlines

    # results statement
    header4 = "Results:"
    result_buffer = summary_buffer - (line_count * 8)
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, result_buffer, header4)
    text_width = c.stringWidth(header4, "Helvetica-Bold", 8)
    c.line(50, result_buffer - 2, 50 + text_width, result_buffer - 2)
    c.setFont("Helvetica", 8)
    result_statement = f"See next page"
    c.drawString(50, result_buffer - 12, result_statement)

    # methods statement
    placeholder_dict = {
        "sample_type": "py_sampletype_ph",
        "na_type": "py_natype_ph",
        "ext_kit": "py_extkit_ph",
        "library_prep": "py_libraryprep_ph",
        "seq_platform": "py_seqplatform_ph",
    }
    
    methods_buffer = result_buffer - (8*6)
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, methods_buffer, header3)
    text_width = c.stringWidth(header3, "Helvetica-Bold", 8)
    c.line(50, methods_buffer - 2, 50 + text_width, methods_buffer - 2)

    methods_buffer = methods_buffer - 15
    for key, value in placeholder_dict.items():
        globals()[key] = in_dict.get(value, value)
        
    methods_statement = f"{sample_type} {na_type} was extracted by {ext_kit}, libraries were generated using the {library_prep}. The library was then sequenced on the {seq_platform}. The sequencing data was analysed using an in-house rapid metagenomics pipeline mirrored from CZ ID (v0.9.1)."
    numlines = draw_wrapped_text(c, methods_statement, methods_buffer, 8)

    # Signatories
    signatories = [['Dr. John Sebastian Eden', 'Dr. Winkie Fong'],
                    ["Post-Doctoral Fellow", "Bioinformatician"],
                    ['Prof Dominic Dwyer', 'Prof Vitali Sintchenko', "Dr. Jen Kok"],
                    ["NSW Health Pathology", "Supervising Pathologist", "Supervising Pathologist"]]
    table = Table(signatories, rowHeights=4*3.5)
    style = TableStyle([
        ('FONTNAME', (0, 0), (2, 0), 'Helvetica-Bold'),
        ('FONTNAME', (0, 2), (2, 2), 'Helvetica-Bold'),
        ('FONTNAME', (0, 1), (2, 1), 'Helvetica'),
        ('FONTNAME', (0, 3), (3, 3), 'Helvetica'),
        ('FONTSIZE', (0, 0), (4, 4), 8),
        ('RIGHTPADDING', (0, 0), (5, -1), 40),
        ('BOTTOMPADDING', (0, 0), (1, 1), 10),
    ])
    table.setStyle(style)
    table.wrapOn(c, 0,0)
    table.drawOn(c, A4[0] / 2 - 180, A4[1] - 700)

    # Disclaimer
    disclaimer_buffer = A4[1] - 750
    c.setFont("Helvetica-Bold", 4)
    disclaimer = f"Disclaimer: The results of this report are confidential and are not to be used or disclosed to any other person or used for any other purpose, unless that use is disclosed or the purpose is expressly authorized in writing by NSW Health Pathology or the named recipient on this report. This test is not accredited by any official certification body. Discretion is advised, as the information presented in the test may not cover all possible medical conditions, and the results may not be exhaustive or up-to-date."
    draw_wrapped_text(c, disclaimer, disclaimer_buffer, 6)
    

    # finish page 1
    c.showPage()


    # start page two for figures
    if in_dict['py_coverageimg_ph']:

        c.setFont("Helvetica", 8)
        figure1 = f"Figure 1: Whole genome coverage map of sequencing reads to {used_reference}"
        c.drawString(60, A4[1] - 60, figure1)
    
        coverage_img = utils.ImageReader(in_dict.get("py_coverageimg_ph", ""))
        c.drawImage(coverage_img, 70, A4[1] / 2 - 50, width=300, height=100, mask="auto")

    # finished
    c.save()

def draw_wrapped_text(c, text, y, font_size):
    # Set the font and font size
    x = 50
    width = 500
    font_name="Helvetica"
    
    c.setFont(font_name, font_size)

    # Create a TextObject for more control
    text_object = c.beginText(x, y)
    text_object.setFont(font_name, font_size)

    # Split the text into lines manually to fit within the specified width
    lines = []
    current_line = ""
    words = text.split()
    
    for word in words:
        text_width = c.stringWidth(current_line + " " + word, font_name, font_size)
        if text_width < width:
            current_line += " " + word
        else:
            lines.append(current_line.strip())
            current_line = word
    
    # Add the last line
    lines.append(current_line.strip())

    # Draw each line
    for line in lines:
        text_object.textLine(line)

    # Draw the text on the canvas
    c.drawText(text_object)

    # Draw each line
    for line in lines:
        text_object.textLine(line)

    # Draw the text on the canvas
    c.drawText(text_object)
    return len(lines)