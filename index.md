---
layout: workshop
venue: "University at Buffalo, SUNY"   # brief name of host site without address (e.g., "Euphoric State University")
address: "University at Buffalo, SUNY, North Campus"     # full street address of workshop (e.g., "Room A, 123 Forth Street, Blimingen, Euphoria")
country: "United States"      # lowercase two-letter ISO country code such as "fr" (see https://en.wikipedia.org/wiki/ISO_3166-1#Current_codes)
language: "en"     # lowercase two-letter ISO language code such as "fr" (see https://en.wikipedia.org/wiki/List_of_ISO_639-1_codes)
latitude: "43.002890"     # decimal latitude of workshop venue (use https://www.latlong.net/)
longitude: "-78.788780"    # decimal longitude of the workshop venue (use https://www.latlong.net)
humandate: "Aug 2-16, 2020"    # human-readable dates for the workshop (e.g., "Feb 17-18, 2020")
humantime: "9:00 am - 4:30 pm"    # human-readable times for the workshop (e.g., "9:00 am - 4:30 pm")
startdate: 2020-08-02      # machine-readable start date for the workshop in YYYY-MM-DD format like 2015-01-01
enddate: 2020-08-16        # machine-readable end date for the workshop in YYYY-MM-DD format like 2015-01-02
instructor: ["Alexey Akimov", "Jeanette Sperhac", "Ivan Infante", 
             "Sergei Tretiak", "Leticia Gonzalez", "Markus Oppel",
             "Sebastian Mai"]  # boxed, comma-separated list of instructors' names as strings, like ["Kay McNulty", "Betty Jennings", "Betty Snyder"]
helper: ["Brendan Smith"]     # boxed, comma-separated list of helpers' names, like ["Marlyn Wescoff", "Fran Bilas", "Ruth Lichterman"]
email: ["alexeyak@buffalo.edu"]    # boxed, comma-separated list of contact email addresses for the host, lead instructor, or whoever else is handling questions, like ["marlyn.wescoff@example.org", "fran.bilas@example.org", "ruth.lichterman@example.org"]
collaborative_notes:         # optional: URL for the workshop collaborative notes, e.g. an Etherpad or Google Docs document (e.g., https://pad.carpentries.org/2015-01-01-euphoria)
googleform: https://forms.gle/bppFcFe1JkCWCu987
carpentry: "sc"
---


{% comment %}
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
{% endcomment %}


> # UPDATE: Due to COVID-19 situation, the Summer School is cancelled for the summer 2020. We plan to run it in summer 2021. Please stay tuned for further updates.
> {: .testimonial}

# Excited States and Nonadiabatic Dynamics Summer School 

## About the Summer School

The summer school aims to educate graduate students, postdocs, and researchers working in a broader field
of nonadiabatic and excited state dynamics as well as in  computational material sciences.
The summer school will provide conceptual and practical hands-on training in a range of methods and cyberinfrastructure tools
for modeling excited states dynamics in abstract models and in atomistic materials.

Students will not only learn about using the tools and the underlying theory, but also how these tools can facilitate their ongoing research
The programming-driven nature of the school will contribute to training the new generation workforce with skills that go beyond
the standard computational chemistry curriculum and are highly desired by both industry and academia.


Keywords and topics:
- nonadiabatic dynamics
- excited states
- quantum dynamics
- quantum-classical methods
- charge transfer
- excitation energy transfer
- trajectory surface hopping
- TD-DFT
- algorithms and methods
- software, programming, Python
- best practices, git, github

 
The school aims to provide training in a range of advanced tools. This year, the focus will be on the following packages:
- Libra (Akimov)
- NEXMD (Tretiak)
- SHARC (Gonzalez, Mai)
- QMflows and QMflows-NAMD (Infante)
- DFTB+
- CP2K
- Quantum Espresso
- ErgoSCF
- OpenMolcas

The school will leverage the VIDIA gateway at University at Buffalo


## Logistics

{% if page.humandate %}
<p id="when">
  <strong>When:</strong>
  {{page.humandate}}.
  {% include workshop_calendar.html %}
</p>
{% endif %}

{% if page.latitude and page.longitude %}
<p id="where">
  <strong>Where:</strong>
  {{page.address}}.
  Get directions with
  <a href="//www.openstreetmap.org/?mlat={{page.latitude}}&mlon={{page.longitude}}&zoom=16">OpenStreetMap</a>
  or
  <a href="//maps.google.com/maps?q={{page.latitude}},{{page.longitude}}">Google Maps</a>.
</p>
{% endif %}

{% comment %}
CONTACT EMAIL ADDRESS
Display the contact email address set in the configuration file.
{% endcomment %}
<p id="contact">
  <strong>Contact</strong>:
  Please email
  {% if page.email %}
  {% for email in page.email %}
  {% if forloop.last and page.email.size > 1 %}
  or
  {% else %}
  {% unless forloop.first %}
  ,
  {% endunless %}
  {% endif %}
  <a href='mailto:{{email}}'>{{email}}</a>
  {% endfor %}
  {% else %}
  to-be-announced
  {% endif %}
  for more information.
</p>

### Schedule

{% include base_path.html %}

The details may vary and some topics may be changed, omitted, or their order may be changed. Please check for updates. 

  <table class="table table-striped">
  <tr>
    <td class="col-md-3">August 2, 2020 (Day 1), <strong>Sunday</strong></td>
    <td class="col-md-7">Arrival. check-in. welcome dinner</td>
    <td class="col-md-2">Schedule of Day 1</td>
  </tr>

  <tr>
    <td class="col-md-3">August 3, 2020 (Day 2), Monday</td>
    <td class="col-md-7">Opening. <a href="01-introduction">Introduction. VIDIA gateway. Setting up.  
    Revision of Python programming. Revision of the best practices</a></td>
    <td class="col-md-2">Schedule of Day 2</td>      
  </tr>

  <tr>
    <td class="col-md-3">August 4, 2020 (Day 3), Tuesday</td>
    <td class="col-md-7"><a href="02-methodology-libra">Methodology prototyping in Python. 
    Introduction of the Libra package. Classical molecular dynamics and convenience functions.</a></td>
    <td class="col-md-2">Schedule of Day 3</td>
  </tr>


  <tr>
    <td class="col-md-3">August 5, 2020 (Day 4), Wednesday</td>
    <td class="col-md-7">Wavepacket dynamics on a grid. Hierarchy of equations of motion (HEOM).</td>
    <td class="col-md-2">Schedule of Day 4</td>
  </tr>

  <tr>
    <td class="col-md-3">August 6, 2020 (Day 5), Thursday</td>
    <td class="col-md-7">Trajectory surface hopping and Ehrenfest methods for model Hamiltonians.</td>
    <td class="col-md-2">Schedule of Day 5</td>
  </tr>

  <tr>
    <td class="col-md-3">August 7, 2020 (Day 6), Friday</td>
    <td class="col-md-7">Nonadiabatic dynamics of atomistics systems: Libra/ErgoSCF, Libra/DFTB+, and Libra/QuantumEspresso.</td>
    <td class="col-md-2">Schedule of Day 6</td>
  </tr>

  <tr>
    <td class="col-md-3">August 8, 2020 (Day 7), Saturday</td>
    <td class="col-md-7">On your own. Projects time</td>
    <td class="col-md-2">Schedule of Day 7</td>
  </tr>

  <tr>
    <td class="col-md-3">August 9, 2020 (Day 8), <strong>Sunday</strong></td>
    <td class="col-md-7">On your own. Projects time</td>
    <td class="col-md-2">Schedule of Day 8</td>
  </tr>


  <tr>
    <td class="col-md-3">August 10, 2020 (Day 9), Monday</td>
    <td class="col-md-7">NA-MD of nanoclusters: Hands on with QMflows-NAMD and CP2K</td>
    <td class="col-md-2">Schedule of Day 9</td>
  </tr>

  <tr>
    <td class="col-md-3">August 11, 2020 (Day 10), Tuesday</td>
    <td class="col-md-7">NA-MD of nanoclusters and condensed matter: Hands on with Libra workflows, DFTB+ and Quantum Espresso</td>
    <td class="col-md-2">Schedule of Day 10</td>
  </tr>

  <tr>
    <td class="col-md-3">August 12, 2020 (Day 11), Wednesday</td>
    <td class="col-md-7">NA-MD of large organic molecules: Hands on with NEXMD</td>
    <td class="col-md-2">Schedule of Day 11</td>
  </tr>

  <tr>
    <td class="col-md-3">August 13, 2020 (Day 12), Thursday</td>
    <td class="col-md-7">Wavepacket dynamics tutorial. NA-MD of small molecules: Hands on with SHARC and OpenMolcas</td>
    <td class="col-md-2">Schedule of Day 12</td>
  </tr>

  <tr>
    <td class="col-md-3">August 14, 2020 (Day 13), Friday</td>
    <td class="col-md-7"> ... Continued </td>
    <td class="col-md-2">Schedule of Day 13</td>
  </tr>

  <tr>
    <td class="col-md-3">August 15, 2020 (Day 14), Saturday</td>
    <td class="col-md-7">Project presentations</td>
    <td class="col-md-2">Schedule of Day 14</td>
  </tr>

  <tr>
    <td class="col-md-3">August 16, 2020 (Day 15), <strong>Sunday</strong></td>
    <td class="col-md-7">Depart</td>
    <td class="col-md-2"></td>
  </tr>
  </table>



## Participation
### How to apply to the school

1. Read this page carefully
2. Prepare your application package (you will need it in the next step)

   2.1. your CV (including graduate GPA)

   2.2. a statement of purpose PDF should describe in no more than 2 pages:

   * why you wish to attend the school
   * your current research project(s)
   * how you think this school will help with your career
   * The summer school will culminate in an elaboration of the participant-designed projects that
     leverage some of the cyberinfrastructure that will be covered in the workshop.
     Please make sure to mention briefly your potential project and how you could leverage any of the
     packages this workshop will teach you.
         
   2.3. request your adviser to submit a letter of recommendation for you to "alexeyak AT buffalo DOT edu", please replace "AT" and "DOT"
   with the corresponding characters

3. Complete the <a href="https://forms.gle/bppFcFe1JkCWCu987" target="_blank" rel="nofollow">**Registration form**</a>


### Important dates
   * School application materials are due by June 1 , 2020
   * Students and Postdocs will be notified of their admission by June 30, 2020


### Who can apply

This summer school is primarily for graduate students working in computational 
modeling of excited states and nonadiabatic dynamics, both in abstract and atomistic
applications/problems. The school aims to help researchers/students working either in 
methodology development for nonadiabatic or quantum-classical dynamics and in 
applied studies of various types of solar energy materials (photovoltaics, photocatalytics, etc.). 

Postdocs and researchers wishing to acquire the practical experience with new simulation
tools and expand their knowledge in the areas of excited states and nonadiabatic dynamics
are also welcomed to participate.




### Selection and support

The applicants will be screened based on several factors, to ensure the adequate level of prior preparation.

Thanks to generous support of the NSF-OAC Cybertraining program, we will cover the travel/lodging expenses of the selected students
up to a reasonable amount.

We can support only a limited number of participants. However, if you are not selected for the support, you can 
still be invited to attend the event with the understanding that you would cover all of your (travel and lodging) expenses.
Please indicate in the registration form if you would like to be considered for self-supported participation.

Out of all applicants, the total of **30 students** will be selected on a competitive basis.  





{% include links.md %}
