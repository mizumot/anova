library(shiny)


shinyUI(pageWithSidebar(


    headerPanel("ANOVA"),


    sidebarPanel(

        p('This web application is completely based on the source code of',
        a("anovakun.", href="http://riseki.php.xdomain.jp/index.php?ANOVA%E5%90%9B", target="_blank"),
        ''),

        br(),

        p(HTML("<hr>")),


        radioButtons("factor", strong("Factors:"),
                list("One-way ANOVA" = "oneway",
                    "Two-way ANOVA" = "twoway"), 'Two-way ANOVA'),
                conditionalPanel(
                        condition = "input.factor == 'oneway'",
                        selectInput("one.design", "Design:",
                        list("Between", "Within"))),
                conditionalPanel(
                        condition = "input.factor == 'twoway'",
                        selectInput("two.design", "Design:",
                        list("Factor1Between_Factor2Between", "Factor1Between_Factor2Within", "Factor1Within_Factor2Within"), 'Factor1Between_Factor2Within')
                ),
        br(),

        p(strong("Levels")),

        numericInput("factor1.level", "Number of levels in Factor 1:", 2),

        conditionalPanel(condition = "input.factor == 'twoway'",
                numericInput("factor2.level", "Number of levels in Factor 2:", 3)
                ),

        br()

    ),


mainPanel(



tabsetPanel(

        tabPanel("Main",

            p('Note: Input values must be separated by tabs. Copy and paste from Excel/Numbers.'),

            p(HTML("<b><div style='background-color:#FADDF2;border:1px solid black;'>Please make sure that your data includes the header (variable names) in the first row.</div></b>")),

            aceEditor("text", value="Method\tPre\tPost\tDelayed\n1\t31\t48\t30\n1\t39\t51\t44\n1\t56\t67\t58\n1\t47\t44\t50\n1\t29\t33\t47\n1\t37\t41\t43\n1\t46\t43\t55\n1\t37\t53\t42\n1\t38\t64\t49\n1\t30\t52\t33\n1\t33\t53\t43\n1\t33\t44\t40\n1\t31\t39\t44\n1\t25\t32\t31\n1\t51\t62\t57\n1\t31\t43\t38\n1\t56\t59\t59\n1\t18\t19\t22\n1\t35\t46\t37\n1\t30\t50\t35\n1\t46\t62\t62\n1\t35\t45\t43\n1\t43\t58\t51\n1\t40\t49\t53\n1\t46\t58\t51\n1\t50\t66\t69\n1\t39\t44\t54\n1\t45\t64\t44\n1\t22\t45\t41\n1\t33\t53\t44\n2\t36\t42\t31\n2\t39\t41\t38\n2\t39\t44\t42\n2\t42\t30\t30\n2\t17\t13\t27\n2\t38\t32\t32\n2\t39\t27\t33\n2\t35\t41\t38\n2\t39\t55\t55\n2\t29\t42\t35\n2\t32\t43\t44\n2\t53\t49\t46\n2\t44\t39\t43\n2\t45\t38\t34\n2\t40\t43\t44\n2\t33\t35\t34\n2\t24\t24\t27\n2\t42\t28\t23\n2\t35\t36\t31\n2\t45\t51\t45\n2\t36\t45\t54\n2\t40\t39\t39\n2\t15\t26\t31\n2\t37\t37\t44\n2\t21\t29\t31\n2\t52\t58\t68\n2\t35\t31\t41\n2\t55\t62\t50\n2\t55\t60\t64\n2\t59\t64\t61",
                mode="r", theme="cobalt"),

            br(),

            h3("Output"),
            verbatimTextOutput("anovakun.out"),

            br(),
            h3("Plot"),
            downloadButton('downloadANOVAPlot1', 'Download the Factor 1 plot as pdf'), br(),

            conditionalPanel(condition = "input.factor == 'twoway'",
                downloadButton('downloadANOVAPlot2', 'Download the Factor 2 plot as pdf')
            ),

            radioButtons("axis", "y-axis",
                list("Default" = "default", "Min-Max" = "min.max"), 'Default'),
            plotOutput("AnovaPlot1", width="80%"),
            plotOutput("AnovaPlot2", width="80%")

            ),



tabPanel("Input Examples",

p('Note: Input values must be separated by tabs. Copy and paste from Excel/Numbers.'),

p(HTML("<b><div style='background-color:#FADDF2;border:1px solid black;'>Please make sure that your data includes the header (variable names) in the first row.</div></b>")),

br(),

p(strong("One-way ANOVA (Between Subjects [3 levels])")),
aceEditor("text1", value="Class\tScore\n1\t76\n1\t54\n1\t62\n1\t46\n1\t53\n1\t64\n1\t42\n1\t96\n1\t87\n1\t92\n1\t89\n1\t63\n1\t39\n1\t30\n1\t85\n1\t31\n1\t60\n1\t63\n1\t94\n1\t29\n1\t47\n1\t53\n1\t33\n1\t56\n1\t46\n1\t97\n1\t54\n1\t43\n1\t79\n2\t56\n2\t70\n2\t88\n2\t42\n2\t64\n2\t40\n2\t62\n2\t86\n2\t55\n2\t45\n2\t33\n2\t77\n2\t59\n2\t82\n2\t80\n2\t95\n2\t89\n2\t35\n2\t92\n2\t88\n2\t53\n2\t43\n2\t47\n2\t94\n2\t64\n2\t78\n2\t29\n2\t47\n2\t65\n3\t57\n3\t39\n3\t84\n3\t63\n3\t84\n3\t54\n3\t86\n3\t50\n3\t93\n3\t30\n3\t63\n3\t66\n3\t46\n3\t35\n3\t98\n3\t86\n3\t85\n3\t73\n3\t98\n3\t99\n3\t39\n3\t50\n3\t87\n3\t67\n3\t57\n3\t80\n3\t92\n3\t41\n3\t82", mode="r", theme="solarized_light"),


br(),
p(strong("One-way ANOVA (Within Subjects [3 levels])")),
aceEditor("text2", value="First\tSecond\tThird\n76\t56\t57\n54\t70\t39\n62\t88\t84\n46\t42\t63\n53\t64\t84\n64\t40\t54\n42\t62\t86\n96\t86\t50\n87\t55\t93\n92\t45\t30\n89\t33\t63\n63\t77\t66\n39\t59\t46\n30\t82\t35\n85\t80\t98\n31\t95\t86\n60\t89\t85\n63\t35\t73\n94\t92\t98\n29\t88\t99\n47\t53\t39\n53\t43\t50\n33\t47\t87\n56\t94\t67\n46\t64\t57\n97\t78\t80\n54\t29\t92\n43\t47\t41\n79\t65\t82", mode="r", theme="solarized_light"),


br(),
p(strong("Two-way ANOVA (Between Subjects [2 levels] & Between Subjects [3 levels])")),
aceEditor("text3", value="Classroom\tClassSize\tScore\n1\t1\t57\n1\t1\t69\n1\t1\t71\n1\t1\t69\n1\t1\t86\n1\t1\t70\n1\t1\t78\n1\t1\t69\n1\t1\t79\n1\t1\t63\n1\t1\t76\n1\t1\t71\n1\t1\t79\n1\t1\t66\n1\t1\t75\n1\t1\t70\n1\t1\t72\n1\t1\t61\n1\t1\t63\n1\t1\t67\n1\t1\t89\n1\t1\t73\n1\t1\t76\n1\t1\t83\n1\t1\t72\n1\t1\t95\n1\t1\t85\n1\t1\t62\n1\t1\t86\n1\t1\t77\n1\t2\t61\n1\t2\t57\n1\t2\t52\n1\t2\t54\n1\t2\t67\n1\t2\t52\n1\t2\t80\n1\t2\t68\n1\t2\t73\n1\t2\t62\n1\t2\t62\n1\t2\t68\n1\t2\t76\n1\t2\t56\n1\t2\t75\n1\t2\t61\n1\t2\t66\n1\t2\t67\n1\t2\t72\n1\t2\t72\n1\t2\t58\n1\t2\t70\n1\t2\t48\n1\t2\t66\n1\t2\t58\n1\t2\t64\n1\t2\t56\n1\t2\t68\n1\t2\t73\n1\t2\t61\n1\t3\t41\n1\t3\t44\n1\t3\t44\n1\t3\t47\n1\t3\t22\n1\t3\t43\n1\t3\t44\n1\t3\t40\n1\t3\t44\n1\t3\t34\n1\t3\t37\n1\t3\t58\n1\t3\t49\n1\t3\t50\n1\t3\t45\n1\t3\t38\n1\t3\t29\n1\t3\t47\n1\t3\t40\n1\t3\t50\n1\t3\t41\n1\t3\t45\n1\t3\t20\n1\t3\t42\n1\t3\t26\n1\t3\t57\n1\t3\t40\n1\t3\t60\n1\t3\t60\n1\t3\t64\n2\t1\t69\n2\t1\t56\n2\t1\t71\n2\t1\t59\n2\t1\t57\n2\t1\t61\n2\t1\t54\n2\t1\t63\n2\t1\t70\n2\t1\t56\n2\t1\t47\n2\t1\t44\n2\t1\t81\n2\t1\t79\n2\t1\t67\n2\t1\t71\n2\t1\t59\n2\t1\t61\n2\t1\t61\n2\t1\t62\n2\t1\t71\n2\t1\t66\n2\t1\t62\n2\t1\t69\n2\t1\t64\n2\t1\t55\n2\t1\t68\n2\t1\t69\n2\t1\t70\n2\t1\t75\n2\t2\t53\n2\t2\t56\n2\t2\t61\n2\t2\t47\n2\t2\t50\n2\t2\t72\n2\t2\t61\n2\t2\t64\n2\t2\t65\n2\t2\t62\n2\t2\t57\n2\t2\t80\n2\t2\t71\n2\t2\t54\n2\t2\t68\n2\t2\t45\n2\t2\t45\n2\t2\t72\n2\t2\t57\n2\t2\t62\n2\t2\t61\n2\t2\t54\n2\t2\t65\n2\t2\t58\n2\t2\t75\n2\t2\t72\n2\t2\t85\n2\t2\t71\n2\t2\t70\n2\t2\t73\n2\t3\t50\n2\t3\t42\n2\t3\t61\n2\t3\t88\n2\t3\t54\n2\t3\t66\n2\t3\t53\n2\t3\t31\n2\t3\t61\n2\t3\t25\n2\t3\t68\n2\t3\t50\n2\t3\t54\n2\t3\t71\n2\t3\t46\n2\t3\t51\n2\t3\t61\n2\t3\t53\n2\t3\t75\n2\t3\t49\n2\t3\t72\n2\t3\t56\n2\t3\t63\n2\t3\t50\n2\t3\t77\n2\t3\t54\n2\t3\t39\n2\t3\t62\n2\t3\t57\n2\t3\t60",mode="r", theme="solarized_light"),


br(),
p(strong("Two-way ANOVA (Between Subjects [2 levels] & Within Subjects [3 levels])")),
aceEditor("text4", value="Method\tPre\tPost\tDelayed\n1\t31\t48\t30\n1\t39\t51\t44\n1\t56\t67\t58\n1\t47\t44\t50\n1\t29\t33\t47\n1\t37\t41\t43\n1\t46\t43\t55\n1\t37\t53\t42\n1\t38\t64\t49\n1\t30\t52\t33\n1\t33\t53\t43\n1\t33\t44\t40\n1\t31\t39\t44\n1\t25\t32\t31\n1\t51\t62\t57\n1\t31\t43\t38\n1\t56\t59\t59\n1\t18\t19\t22\n1\t35\t46\t37\n1\t30\t50\t35\n1\t46\t62\t62\n1\t35\t45\t43\n1\t43\t58\t51\n1\t40\t49\t53\n1\t46\t58\t51\n1\t50\t66\t69\n1\t39\t44\t54\n1\t45\t64\t44\n1\t22\t45\t41\n1\t33\t53\t44\n2\t36\t42\t31\n2\t39\t41\t38\n2\t39\t44\t42\n2\t42\t30\t30\n2\t17\t13\t27\n2\t38\t32\t32\n2\t39\t27\t33\n2\t35\t41\t38\n2\t39\t55\t55\n2\t29\t42\t35\n2\t32\t43\t44\n2\t53\t49\t46\n2\t44\t39\t43\n2\t45\t38\t34\n2\t40\t43\t44\n2\t33\t35\t34\n2\t24\t24\t27\n2\t42\t28\t23\n2\t35\t36\t31\n2\t45\t51\t45\n2\t36\t45\t54\n2\t40\t39\t39\n2\t15\t26\t31\n2\t37\t37\t44\n2\t21\t29\t31\n2\t52\t58\t68\n2\t35\t31\t41\n2\t55\t62\t50\n2\t55\t60\t64\n2\t59\t64\t61", mode="r", theme="solarized_light"),


br(),
p(strong("Two-way ANOVA (Witin Subjects[2 levels] * Within Subjects [3 levels])")),
aceEditor("text5", value="TypeA.Pre\tTypeA.Post\tTypeA.Delayed\tTypeB.Pre\tTypeB.Post\tTypeB.Delayed\n44\t120\t153\t51\t100\t110\n61\t119\t148\t62\t109\t117\n67\t157\t167\t56\t134\t139\n60\t153\t175\t57\t140\t161\n61\t139\t162\t59\t126\t137\n55\t130\t161\t67\t120\t155\n73\t135\t172\t61\t122\t162\n55\t122\t144\t73\t115\t145\n62\t144\t165\t62\t110\t162\n60\t155\t172\t52\t103\t159",
mode="r", theme="solarized_light"),

br()

),






        tabPanel("About",

            strong('Note'),
            p('This web application is developed with',
            a("Shiny.", href="http://www.rstudio.com/shiny/", target="_blank"),
            ''),

            br(),

            strong('List of Packages Used'), br(),
            code('library(shiny)'),br(),
            code('library(shinyAce)'),br(),
            code('library(sciplot)'),br(),

            br(),

            strong('Code'),
            p('Examples are based on',
            a('"The handbook of Research in Foreign Language Learning and Teaching" (Takeuchi & Mizumoto, 2012).', href='http://mizumot.com/handbook/', target="_blank")),

            p('Source code for this application is based on',
            a('anovakun.', href="http://riseki.php.xdomain.jp/index.php?ANOVA%E5%90%9B", target="_blank")),

            p('The code for this web application is available at',
            a('GitHub.', href='https://github.com/mizumot/anova', target="_blank")),

            p('If you want to run this code on your computer (in a local R session), run the code below:',
            br(),
            code('library(shiny)'),br(),
            code('runGitHub("anova","mizumot")')
            ),

            br(),

            strong('Recommended'),
            p('To learn more about R, I suggest this excellent and free e-book (pdf),',
            a("A Guide to Doing Statistics in Second Language Research Using R,", href="http://cw.routledge.com/textbooks/9780805861853/guide-to-R.asp", target="_blank"),
            'written by Dr. Jenifer Larson-Hall.'),

            p('Also, if you are a cool Mac user and want to use R with GUI,',
            a("MacR", href="http://www.urano-ken.com/blog/2013/02/25/installing-and-using-macr/", target="_blank"),
            'is defenitely the way to go!'),

            br(),

            strong('Author'),
            p(a("Atsushi MIZUMOTO,", href="http://mizumot.com", target="_blank"),' Ph.D.',br(),
            'Associate Professor of Applied Linguistics',br(),
            'Faculty of Foreign Language Studies /',br(),
            'Graduate School of Foreign Language Education and Research,',br(),
            'Kansai University, Osaka, Japan'),

            br(),

            a(img(src="http://i.creativecommons.org/p/mark/1.0/80x15.png"), target="_blank", href="http://creativecommons.org/publicdomain/mark/1.0/"),

            p(br())

)
)
)
))
