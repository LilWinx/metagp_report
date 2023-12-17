from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib import utils
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Table, TableStyle, Paragraph
import numpy as np

def pdf_template(outdir, in_dict, used_reference, img_type):
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
    # coverage figure if exists
    figure_count = 1
    coverage = False
    shift_down = 0
    c.setFont("Helvetica", 8)
    if in_dict['py_coverageimg_ph']:
        coverage = True
        figure1 = f"Figure {figure_count}: Whole genome coverage map of sequencing reads to {used_reference}"
        figure_count += 1
        c.drawString(60, A4[1] - 60, figure1)
    
        coverage_img = utils.ImageReader(in_dict.get("py_coverageimg_ph", ""))
        c.drawImage(coverage_img, 60, A4[1] - 220, width=500, height=150, mask="auto")
        shift_down = 200
    
    # hbar image.
    if img_type == "hbar":
        figure2 = f"Figure {figure_count}: Relative abundance of the Top 10 species following filtering."
        figure_count += 1
        
        c.drawString(60, A4[1] - 60 - shift_down, figure2)
        hbar_img = utils.ImageReader(in_dict.get("py_hbar_ph", ""))
        c.drawImage(hbar_img, 60, A4[1] - 180 - shift_down, width=500, height=100, mask="auto")
    elif img_type == "donut":
        figure2 = f"Figure {figure_count}: Relative abundance of the Top 10 species following filtering from both DNA and RNA mNGS."
        figure_count += 1
        
        c.drawString(60, A4[1] - 60 - shift_down, figure2)
        hbar_img = utils.ImageReader(in_dict.get("py_hbar_ph", ""))
        c.drawImage(hbar_img, 60, A4[1] - 520 - shift_down, width=400, height=450, mask="auto")
    else:
        pass

    # end page 2
    c.showPage()
    
    # Top10 species table
    shift_down = 0
    c.setFont("Helvetica-Bold", 8)
    c.drawString(50, A4[1] - 60, f"Top 10 Species")
    
    top10 = [
        'py_tspecies1_ph',
        'py_tspecies1status_ph',
        'py_tspecies1zscore_ph',
        'py_tspecies2_ph',
        'py_tspecies2status_ph',
        'py_tspecies3zscore_ph',
        'py_tspecies3_ph',
        'py_tspecies3status_ph',
        'py_tspecies3zscore_ph',
        'py_tspecies4_ph',
        'py_tspecies4status_ph',
        'py_tspecies4zscore_ph',
        'py_tspecies5_ph',
        'py_tspecies5status_ph',
        'py_tspecies5zscore_ph',
        'py_tspecies6_ph',
        'py_tspecies6status_ph',
        'py_tspecies6zscore_ph',
        'py_tspecies7_ph',
        'py_tspecies7status_ph',
        'py_tspecies7zscore_ph',
        'py_tspecies8_ph',
        'py_tspecies8status_ph',
        'py_tspecies8zscore_ph',
        'py_tspecies9_ph',
        'py_tspecies9status_ph',
        'py_tspecies9zscore_ph',
        'py_tspecies10_ph',
        'py_tspecies10status_ph',
        'py_tspecies10zscore_ph',
    ]

    for text in range(len(top10)):

        for key, value in in_dict.items():
            # Check if the placeholder contains "zscore"
            value_str = str(value)
            value_str = value_str.replace('❗', ' !!')
            if "zscore" in key and key in top10[text]:
                top10[text] = top10[text].replace(key, f'Z-score: {value_str}')
            else:
                
                top10[text] = top10[text].replace(key, value_str)
    
    
    table_buffer = A4[1] - 75
    for i, text in enumerate(top10):
        if i % 3 == 0 or i == 0:
            c.setFont("Helvetica-Bold", 8)
            c.drawString(50, table_buffer, text)
            numlines = 1
        else:
            numlines = draw_wrapped_text(c, text, table_buffer, 8)
        table_buffer -= 12 * numlines
    """
    note = f"Note: A Z-Score of 100 is considered statistically significant, the closer the Z-Score is to 1, the species is considered insignificant, or was found within the negative control and considered background species."
    draw_wrapped_text(c, note, table_buffer - 5, 6)
    """

    c.showPage()

    
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